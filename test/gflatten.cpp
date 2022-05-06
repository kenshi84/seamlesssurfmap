#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>

#include <vector>

#include "util.h"
#include "mosek_util.h"

struct SeamVertex {
  Vertex left, right;
  Vector2d alpha;
  Variable::t alphaVar;   // MOSEK variable
};

using SeamVertices = std::vector<SeamVertex>;

struct SeamInfo {
  SeamVertices vertices;

  Vector2d similarityLeft;              // (a, b) treated as a 2x2 matrix [a, -b; b, a], as in Equation (2)
  Vector2d similarityRight;

  // MOSEK variables
  Variable::t similarityLeftVar;
  Variable::t similarityRightVar;
};

// Local per-triangle computation of map differential & best fit rotation
void gflatten_local(const MatrixXd& V, const MatrixXi& F, const MatrixXd& V_uv, std::vector<util::MapDifferential>& mapDiff, std::vector<Matrix2d>& R, std::vector<double>& c) {
  const int nF = F.rows();

  mapDiff.resize(nF);
  R.resize(nF);
  c.resize(nF);

  for (int i = 0; i < nF; ++i) {
    mapDiff[i] = util::compute_map_differential(V, F, V_uv, i);

    Matrix2d Y;
    util::polar_decomposition(mapDiff[i].A, R[i], Y);

    double detA = mapDiff[i].A.determinant();
    c[i] = 1.0 / std::sqrt(std::abs(detA));

    c[i] = std::max<double>(std::min<double>(c[i], 100), 0.01);
  }
}

// Global optimization via SOCP
void gflatten_global(const MatrixXd& V, const MatrixXi& F, const std::vector<util::MapDifferential>& mapDiff, const std::vector<Matrix2d>& R, const std::vector<double>& c, int iter, const std::vector<Vertex>& landmarkRing, std::vector<SeamInfo>& seamInfos, MatrixXd& V_uv) {
  MosekModel::t M = new MosekModel("gflatten");
  auto _M = finally([&]() { M->dispose(); });

  const int nV = V.rows();
  const int nF = F.rows();
  const int nS = (int)seamInfos.size();

  ASSERT(V.cols() == 3);
  ASSERT(F.cols() == 3);
  ASSERT(mapDiff.size() == nF);
  ASSERT(R.size() == nF);
  ASSERT(c.size() == nF);

  V_uv.resize(nV, 2);

  const double sqrt2 = std::sqrt(2.0);

  // Main variable (UV coordinate) per vertex
  Variable::t uv = M->variable("uv", 2 * nV, Domain::unbounded());

  auto get_variable_uv_vector = [&](int vi) {
    ASSERT(vi < nV);
    return Expr::vstack(uv->index(vi), uv->index(nV + vi));
  };

  // Auxiliary variable per triangle
  Variable::t rB = M->variable("rB", nF, Domain::unbounded());
  Variable::t rC = M->variable("rC", nF, Domain::unbounded());
  Variable::t u  = M->variable("u" , nF, Domain::unbounded());
  Variable::t s  = M->variable("s" , nF, Domain::unbounded());
  Variable::t S  = M->variable("S" , nF, Domain::unbounded());
  Variable::t U  = M->variable("U" , nF, Domain::unbounded());

  std::vector<double> sqrt_area(nF);

  SPDLOG_INFO("Building model...");

  for (int i = 0; i < nF; ++i) {
    auto rB_i = rB->index(i);
    auto rC_i = rC->index(i);
    auto u_i  = u ->index(i);
    auto s_i  = s ->index(i);
    auto S_i  = S ->index(i);
    auto U_i  = U ->index(i);

    // Skip degenerate or inverted triangles
    if (mapDiff[i].Q.determinant() <= 0) {
      SPDLOG_WARN("Face {} is degenerate or inverted! Skipping", i);
      M->constraint("degenerate_" + std::to_string(i),
        U_i,
        Domain::equalsTo(0.0));
      continue;
    }

    auto Pinv_MOSEK = mosek_util::from_eigen_to_mosek(mapDiff[i].Pinv);

    sqrt_area[i] = std::sqrt(0.5 * mapDiff[i].P.determinant());

    // Variable 2D vectors representing the triangle corners
    int v0 = F(i, 0);
    int v1 = F(i, 1);
    int v2 = F(i, 2);

    auto q0 = get_variable_uv_vector(v0);
    auto q1 = get_variable_uv_vector(v1);
    auto q2 = get_variable_uv_vector(v2);

    auto q01 = Expr::sub(q1, q0);
    auto q02 = Expr::sub(q2, q0);

    auto Q = Expr::hstack(q01, q02);

    auto A = Expr::mul(Q, Pinv_MOSEK);

    // Decompose A = B + C where B = 0.5 (A - A^T + tr(A)I)
    auto At = Expr::transpose(A);
    auto trA = mosek_util::trace(A);
    auto trA_I = Expr::mul(trA, MosekMatrix::eye(2));
    auto B = Expr::mul(0.5, Expr::sub(Expr::add(A, trA_I), At));
    auto C = Expr::sub(A, B);

    // Equation (10a) --------------+
    //    |B| + |C| <= sqrt(2) r    |
    //------------------------------+
    // |B| <= rB
    M->constraint("qc_10a_B_" + std::to_string(i),
      Expr::vstack(rB_i, Expr::flatten(B)),
      Domain::inQCone());

    // |C| <= rC
    M->constraint("qc_10a_C_" + std::to_string(i),
      Expr::vstack(rC_i, Expr::flatten(C)),
      Domain::inQCone());

    // r = (rB + rC) / sqrt(2)
    auto r = Expr::mul(1.0 / sqrt2, Expr::add(rB_i, rC_i));

    // Equation (10b) --------------------+
    //    |C| <= tr(R^T B) - sqrt(2) u    |
    //------------------------------------+
    auto R_MOSEK = MosekMatrix::dense(mosek_util::from_eigen_to_mosek(R[i]));
    auto RtB = Expr::mul(R_MOSEK->transpose(), B);
    auto tr_RtB = mosek_util::trace(RtB);

    M->constraint("qc_10b_" + std::to_string(i),
      Expr::vstack(Expr::sub(tr_RtB, Expr::mul(sqrt2, u_i)), Expr::flatten(C)),
      Domain::inQCone());

    // Equation (10c) --------------------+
    //    sqrt((u - s)^2 + 4) <= u + s    |
    //------------------------------------+
    M->constraint("qc_10c_" + std::to_string(i),
      Expr::vstack(Expr::add(u_i, s_i), Expr::sub(u_i, s_i), 2.0),
      Domain::inQCone());

    // Equation (10d) ------------------------+
    //    sqrt((s / c)^2 + (c * r)^2) <= S    |
    //----------------------------------------+
    M->constraint("qc_10d_" + std::to_string(i),
      Expr::vstack(S_i, Expr::mul(1.0 / c[i], s_i), Expr::mul(c[i], r)),
      Domain::inQCone());

    // Equation (10e) --------------------------+
    //    sqrt(S^2 + (U - 1/4)^2) <= U + 1/4    |
    //------------------------------------------+
    M->constraint("qc_10e_" + std::to_string(i),
      Expr::vstack(Expr::add(U_i, 0.25), S_i, Expr::sub(U_i, 0.25)),
      Domain::inQCone());
  }

  // Final SOCP constraint ---------------------+
  //    y >= a_1 * U_1^2 + ... + a_n * U_n^2    |
  //--------------------------------------------+
  Variable::t y = M->variable("y", Domain::unbounded());
  M->constraint("qc_y",
    Expr::vstack(0.5, y, Expr::mulElm(U, new_array_ptr(sqrt_area))),
    Domain::inRotatedQCone());

  // Positional constraints for two landmarks to eliminate remaining DoF (translation & rotation)
  int v0 = landmarkRing[0].getIndex();
  int v1 = landmarkRing[landmarkRing.size() / 2].getIndex();
  SPDLOG_INFO("Fixed vertices: {}, {}", v0, v1);

  Vector2d prev_q0 = V_uv.row(v0);
  Vector2d prev_q1 = V_uv.row(v1);

  Vector2d dir = prev_q1 - prev_q0;
  dir = { -dir[1], dir[0] };
  auto q0 = get_variable_uv_vector(v0);
  auto q1 = get_variable_uv_vector(v1);

  M->constraint("fix_uv_0",
    q0,
    Domain::equalsTo(0.0));

  M->constraint("fix_uv_1",
    Expr::dot(Expr::sub(q1, q0), mosek_util::from_eigen_to_mosek(dir)),
    Domain::equalsTo(0.0));

  if (iter % 2) {
    //-------------------------------------------------------------------------------+
    // Odd iteration: fix similarity transformations, optimize seam vertex positions |
    //-------------------------------------------------------------------------------+
    for (int i = 0; i < nS; ++i) {
      const std::string seam_i_str = "seam_" + std::to_string(i);

      auto T_left = mosek_util::get_similarity_matrix(seamInfos[i].similarityLeft);
      auto T_right = mosek_util::get_similarity_matrix(seamInfos[i].similarityRight);

      auto alpha0 = seamInfos[i].vertices[0].alphaVar = M->variable(seam_i_str + "_alpha_0", 2, Domain::unbounded());

      auto q0_left = get_variable_uv_vector(seamInfos[i].vertices[0].left.getIndex());
      auto q0_right = get_variable_uv_vector(seamInfos[i].vertices[0].right.getIndex());

      for (size_t j = 1; j < seamInfos[i].vertices.size(); ++j) {
        auto alphaj = seamInfos[i].vertices[j].alphaVar = M->variable(seam_i_str + "_alpha_" + std::to_string(j), 2, Domain::unbounded());

        auto qj_left = get_variable_uv_vector(seamInfos[i].vertices[j].left.getIndex());
        auto qj_right = get_variable_uv_vector(seamInfos[i].vertices[j].right.getIndex());

        M->constraint(seam_i_str + "_left_" + std::to_string(j),
          Expr::sub(Expr::sub(qj_left, q0_left), Expr::mul(T_left, Expr::sub(alphaj, alpha0))),
          Domain::equalsTo(0.0));

        M->constraint(seam_i_str + "_right_" + std::to_string(j),
          Expr::sub(Expr::sub(qj_right, q0_right), Expr::mul(T_right, Expr::sub(alphaj, alpha0))),
          Domain::equalsTo(0.0));
      }

      // Fix one endpoint of auxiliary seam vertex to (0, 0)
      M->constraint(seam_i_str + "_fix",
        alpha0,
        Domain::equalsTo(0.0));
    }

  } else {
    //------------------------------------------------------------------------------------------+
    // Even iteration: fix auxiliary seam vertex positions, optimize similarity transformations |
    //------------------------------------------------------------------------------------------+
    for (int i = 0; i < nS; ++i) {
      const std::string seam_i_str = "seam_" + std::to_string(i);

      seamInfos[i].similarityLeftVar = M->variable(seam_i_str + "_similarityLeft", 2, Domain::unbounded());
      seamInfos[i].similarityRightVar = M->variable(seam_i_str + "_similarityRight", 2, Domain::unbounded());

      auto T_left = mosek_util::get_variable_similarity_matrix(seamInfos[i].similarityLeftVar);
      auto T_right = mosek_util::get_variable_similarity_matrix(seamInfos[i].similarityRightVar);

      Vector2d alpha0 = seamInfos[i].vertices[0].alpha;

      auto q0_left = get_variable_uv_vector(seamInfos[i].vertices[0].left.getIndex());
      auto q0_right = get_variable_uv_vector(seamInfos[i].vertices[0].right.getIndex());

      for (size_t j = 1; j < seamInfos[i].vertices.size(); ++j) {
        Vector2d alphaj = seamInfos[i].vertices[j].alpha;

        auto alpha0j = mosek_util::from_eigen_to_mosek(Vector2d(alphaj - alpha0));

        auto qj_left = get_variable_uv_vector(seamInfos[i].vertices[j].left.getIndex());
        auto qj_right = get_variable_uv_vector(seamInfos[i].vertices[j].right.getIndex());

        M->constraint(seam_i_str + "_left_" + std::to_string(j),
          Expr::sub(Expr::sub(qj_left, q0_left), Expr::mul(T_left, alpha0j)),
          Domain::equalsTo(0.0));

        M->constraint(seam_i_str + "_right_" + std::to_string(j),
          Expr::sub(Expr::sub(qj_right, q0_right), Expr::mul(T_right, alpha0j)),
          Domain::equalsTo(0.0));
      }
    }
  }

  // Set the objective function to minimize y
  M->objective("obj", ObjectiveSense::Minimize, Expr::mul(0.25, y));

  SPDLOG_INFO("Solving...");

  // Solve
  M->solve();
  SPDLOG_INFO("Primary objective value: " ANSI_BOLD ANSI_CYAN "{}" ANSI_RESET, M->primalObjValue());

  // Retrieve solution values
  ndarray<double, 1> uv_level = *(uv->level());
  for (int i = 0; i < nV; ++i)
    V_uv.row(i) << uv_level[i], uv_level[nV + i];

  // Retrieve auxiliary solution values
  if (iter % 2) {
    // Odd iteration: copy seam vertex positions
    for (int i = 0; i < nS; ++i) {
      for (size_t j = 0; j < seamInfos[i].vertices.size(); ++j) {
        ndarray<double, 1> alphaj_level = *(seamInfos[i].vertices[j].alphaVar->level());
        seamInfos[i].vertices[j].alpha = kt84::vector_cast<2, Vector2d>(alphaj_level);
      }
    }

  } else {
    // Even iteration: copy similarity transformations
    for (int i = 0; i < nS; ++i) {
      ndarray<double, 1> similarityLeft_level = *(seamInfos[i].similarityLeftVar->level());
      ndarray<double, 1> similarityRight_level = *(seamInfos[i].similarityRightVar->level());
      seamInfos[i].similarityLeft = kt84::vector_cast<2, Vector2d>(similarityLeft_level);
      seamInfos[i].similarityRight = kt84::vector_cast<2, Vector2d>(similarityRight_level);
    }
  }
}

int main(int argc, char *argv[]) {
  // Load a manifold surface mesh from file
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;

#if 0
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/A-disk.obj");
  const std::string landmarkGroups_str = "5105,20247;5433;5805,20013,20105;6722,20199;14143,20059;14839,20151;15941;17597";
#endif

#if 0
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/B-disk.obj");
  const std::string landmarkGroups_str = "8713,9202;4828;6294,9048,9108;2669,9170;3288,9075;2977,9138;7543;1643";
#endif

#if 0
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/vase/A-disk.obj");
  const std::string landmarkGroups_str = "111,2565;1111;1224;1526,2532,2547;2431,2618,2641#1268,2594;2518,2619";
#endif

#if 0
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/vase/B-disk.obj");
  const std::string landmarkGroups_str = "1899,2582;9;2165;409,2552,2566;405,2630,2657#1014,2608;2538,2631";
#endif

#if 0
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pretzel-tripletorus/A-disk.obj");
  const std::string landmarkGroups_str = "9,2589;405,2628;409,2552,2573;1899,2654;2165#1014,2606;2538,2629";
#endif

#if 0
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pretzel-tripletorus/A-disk.obj");
  const std::string landmarkGroups_str = "952,2091,2166;954,2041,2127,2150;1507,2077,2107#109,2065;2026,2152;811,2078;2005,2129;1343,2093;1984,2108";
#endif

#if 1
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/hands/A-disk.obj");
  const std::string landmarkGroups_str = "85;407;2703;3192;3803;3903;4649,11544;6538;6722,11421,11607;7058,11488,11502,11578;8012,11442,11675,11717;8349;8376,11462,11639;9037,11522;9491,11432;10759,11471;10767,11452";
#endif

  // For interoperability with libigl
  MatrixXd V = util::get_vertex_position_matrix(*geometry);
  MatrixXi F = mesh->getFaceVertexMatrix<int>();

  // Unknown vertex positions in the parameter space
  MatrixXd V_uv;

  ASSERT(mesh->genus() == 0);
  ASSERT(mesh->nConnectedComponents() == 1);
  ASSERT(mesh->nBoundaryLoops() == 1);

  BoundaryLoop bnd = mesh->boundaryLoop(0);

  const int nV = mesh->nVertices();

  // Split the string into two at '#'
  size_t sharpPos = landmarkGroups_str.find_first_of('#');
  if (sharpPos != std::string::npos) {
    ASSERT(landmarkGroups_str.find_first_of('#', sharpPos + 1) == std::string::npos);
  }
  const std::string normalLandmarkGroups_str = landmarkGroups_str.substr(0, sharpPos);
  const std::string handleLandmarkGroups_str = landmarkGroups_str.substr(sharpPos + 1);

  std::set<Vertex> allLandmarks;

  // Parse landmark groups from string
  auto parse_landmark_groups = [&](const std::string& str) {
    std::vector<std::set<Vertex>> landmarkGroups;

    for (std::string landmarkGroup_str : util::tokenize(str, ';')) {
      std::set<Vertex> landmarkGroup;

      for (std::string landmark_str : util::tokenize(landmarkGroup_str, ',')) {
        int landmark_vid = std::atoi(landmark_str.c_str());
        ASSERT(landmark_vid < mesh->nVertices());

        Vertex v = mesh->vertex(landmark_vid);

        landmarkGroup.insert(v);
        allLandmarks.insert(v);
      }

      landmarkGroups.push_back(landmarkGroup);
    }

    return landmarkGroups;
  };

  std::vector<std::set<Vertex>> normalLandmarkGroups = parse_landmark_groups(normalLandmarkGroups_str);
  std::vector<std::set<Vertex>> handleLandmarkGroups = parse_landmark_groups(handleLandmarkGroups_str);


  // The number of seams is half the number of all landmarks
  const int nS2 = allLandmarks.size();
  const int nS = nS2 / 2;
  ASSERT(nS2 % 2 == 0);

  // Walk along the boundary to extract boundary sides
  std::vector<std::vector<Vertex>> boundarySides(nS2);
  {
    // Starting landmark vertex
    Vertex vStart = *normalLandmarkGroups[0].begin();
    ASSERT(vStart.isBoundary());

    Halfedge heCurr = vStart.halfedge().twin().next();
    ASSERT(!heCurr.isInterior());
    ASSERT(heCurr.tailVertex() == vStart);

    for (size_t i = 0; i < nS2; ++i) {
      SPDLOG_INFO("Landmark vertex along the disk boundary: {}", heCurr.tailVertex());
      do {
        boundarySides[i].push_back(heCurr.tailVertex());
        heCurr = heCurr.next();
      } while (allLandmarks.count(heCurr.tailVertex()) == 0);
    }

    ASSERT(heCurr.tailVertex() == vStart);
  }

  // Initialize V_uv via Tutte embedding with regular polygonal boundary
  {
    std::vector<std::vector<Vector2d>> boundarySidesPos(nS2);

    for (size_t i = 0; i < nS2; ++i) {
      double theta0 = -(2.0 * util::PI() * i) / nS2;
      double theta1 = -(2.0 * util::PI() * (i + 1)) / nS2;

      Vector2d pos0 = { std::cos(theta0), std::sin(theta0) };
      Vector2d pos1 = { std::cos(theta1), std::sin(theta1) };

      boundarySidesPos[i].resize(boundarySides[i].size());

      for (size_t j = 0; j < boundarySides[i].size(); ++j) {
        double a = j / (double)boundarySides[i].size();
        boundarySidesPos[i][j] = (1.0 - a) * pos0 + a * pos1;
      }
    }

    // Flatten these two arrays
    VectorXi bnd2;
    MatrixXd bnd_uv;

    bnd2.resize(bnd.degree(), 1);
    bnd_uv.resize(bnd.degree(), 2);

    int index = 0;
    for (size_t i = 0; i < nS2; ++i) {
      for (size_t j = 0; j < boundarySides[i].size(); ++j) {
        bnd2(index) = boundarySides[i][j].getIndex();
        bnd_uv.row(index) << boundarySidesPos[i][j].transpose();
        ++index;
      }
    }

    // Harmonic parameterization
    igl::harmonic(F, bnd2, bnd_uv, 1, V_uv);

    util::write_flattened_mesh(V_uv, F, "gflatten-init.obj");
  }

  SPDLOG_INFO("Initial state conformal energy: " ANSI_BOLD ANSI_YELLOW "{}" ANSI_RESET, util::compute_conformal_energy(V, F, V_uv));

  // Assign normal landmark group ID
  VertexData<int> vertexLandmarkID(*mesh, -1);
  for (int i = 0; i < (int)normalLandmarkGroups.size(); ++i) {
    for (auto v : normalLandmarkGroups[i])
      vertexLandmarkID[v] = i;
  }

  // Assign handle landmark group ID (offset by nV)
  for (int i = 0; i < (int)handleLandmarkGroups.size(); ++i) {
    for (auto v : handleLandmarkGroups[i])
      vertexLandmarkID[v] = nV + i;
  }

  // Append front of the next boundary vertex array to the end of current seam vertex array, for convenience
  std::vector<std::vector<Vertex>> seamVertexArrays = boundarySides;
  std::vector<Vertex> landmarkRing(nS2);
  for (int i = 0; i < nS2; ++i) {
    seamVertexArrays[i].push_back(seamVertexArrays[(i + 1) % nS2].front());
    landmarkRing[i] = seamVertexArrays[i].front();
  }

  // Find matching pair in seamVertexArrays, 
  std::vector<SeamInfo> seamInfos;
  seamInfos.reserve(nS);
  while (!seamVertexArrays.empty()) {
    for (auto i = seamVertexArrays.begin(); ; ++i) {
      int id0 = vertexLandmarkID[i->front()];
      int id1 = vertexLandmarkID[i->back()];
      ASSERT(id0 > -1);
      ASSERT(id1 > -1);

      auto addSeamInfo = [&](const std::vector<Vertex>& verticesLeft, const std::vector<Vertex>& verticesRight) {
        seamInfos.resize(seamInfos.size() + 1);
        SeamInfo& seamInfo = seamInfos.back();

        size_t nV = verticesLeft.size();
        ASSERT(nV == verticesRight.size());

        seamInfo.vertices.resize(nV);

        for (size_t j = 0; j < nV; ++j) {
          seamInfo.vertices[j].left = verticesLeft[j];
          seamInfo.vertices[j].right = verticesRight[j];
          seamInfo.vertices[j].alpha = { (double)j, 0.0 };
        }
      };

      if (id0 == id1 && (id0 - nV) % 2 == 0) {
        // Seam runs around a handle, starting from and ending at the same landmark group
        std::vector<Vertex> verticesRight = *i;
        seamVertexArrays.erase(i);

        // Look for the corresponding seam on the left
        i = seamVertexArrays.begin();
        for (; ; ++i) {
          if (vertexLandmarkID[i->front()] == id0 + 1 && vertexLandmarkID[i->back()] == id0 + 1) break;
        }
        ASSERT(i != seamVertexArrays.end());

        std::vector<Vertex> verticesLeft = *i;
        seamVertexArrays.erase(i);

        // Orient the left side the same as the right side
        std::reverse(verticesLeft.begin(), verticesLeft.end());

        addSeamInfo(verticesLeft, verticesRight);

        // Move on to finding the next pair
        break;

      } else if (id0 < id1) {
        // Seam on the right found
        std::vector<Vertex> verticesRight = *i;
        seamVertexArrays.erase(i);

        // Look for the corresponding seam on the left
        i = seamVertexArrays.begin();
        for (; ; ++i) {
          if (vertexLandmarkID[i->front()] == id1 && vertexLandmarkID[i->back()] == id0) break;
        }
        ASSERT(i != seamVertexArrays.end());

        std::vector<Vertex> verticesLeft = *i;
        seamVertexArrays.erase(i);

        // Orient the left side the same as the right side
        std::reverse(verticesLeft.begin(), verticesLeft.end());

        addSeamInfo(verticesLeft, verticesRight);

        // Move on to finding the next pair
        break;
      }
    }
  }

  // Iterative optimization
  double prev_conformal_energy;
  kt84::MinSelector<int> bestIter;
  for (int iter = 0; ; ++iter) {
    SPDLOG_INFO("Iteration {}...", iter);

    std::vector<util::MapDifferential> mapDiff;
    std::vector<Matrix2d> R;
    std::vector<double> c;
    SPDLOG_INFO("Local step...");
    gflatten_local(V, F, V_uv, mapDiff, R, c);

    if (iter == 0) c = std::vector<double>(F.rows(), 1.0);     // Initial values should be 1

    SPDLOG_INFO("Global step...");
    gflatten_global(V, F, mapDiff, R, c, iter, landmarkRing, seamInfos, V_uv);

    // Normalize total area
    double total_area = 0;
    for (int j = 0; j < F.rows(); ++j) {
      Vector2d q0 = V_uv.row(F(j, 0));
      Vector2d q1 = V_uv.row(F(j, 1));
      Vector2d q2 = V_uv.row(F(j, 2));
      Matrix2d Q;
      Q << q1 - q0, q2 - q0;
      total_area += 0.5 * Q.determinant();
    }
    V_uv *= 1.0 / std::sqrt(total_area);

    SPDLOG_INFO("Optimized conformal energy: " ANSI_BOLD ANSI_YELLOW "{}" ANSI_RESET, (prev_conformal_energy = util::compute_conformal_energy(V, F, V_uv)));

    util::write_flattened_mesh(V_uv, F, "gflatten-" + std::to_string(iter) + ".obj");

    bestIter.update(prev_conformal_energy, iter);
    SPDLOG_INFO("Best iteration found so far: " ANSI_BOLD ANSI_MAGENTA "{}" ANSI_RESET, bestIter.value);
  }
}
