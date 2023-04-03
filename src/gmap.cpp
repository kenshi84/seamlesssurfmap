#include <igl/doublearea.h>
#include <igl/harmonic.h>

#include <vector>
#include <csignal>

#include "gmap.h"
#include "util.h"
#include "mosek_util.h"
#include "untangle.h"

namespace {

using Vertices = std::vector<Vertex>;

struct VerticesPair : public Pair<Vertices> {
  VertexPair front() const { return { A.front(), B.front() }; }     // Returns copy, not reference
  VertexPair back() const { return { A.back(), B.back() }; }
  void push_back(const VertexPair& p) {
    A.push_back(p.A);
    B.push_back(p.B);
  }
  void reverse() {
    std::reverse(A.begin(), A.end());
    std::reverse(B.begin(), B.end());
  }
  Pair<size_t> size() const { return { A.size(), B.size() }; }
  void resize(const Pair<size_t>& p) {
    A.resize(p.A);
    B.resize(p.B);
  }
};

struct SeamVertex {
  Vertex left, right;
  Vector2d alpha;
  Variable::t alphaVar;   // MOSEK variable
};

using SeamVertexPair = Pair<SeamVertex>;

using SeamVertices = std::vector<SeamVertex>;

struct SeamVerticesPair : public Pair<SeamVertices> {
  SeamVertexPair front() const { return { A.front(), B.front() }; }
  SeamVertexPair back() const { return { A.back(), B.back() }; }
  void resize(const Pair<size_t>& p) {
    A.resize(p.A);
    B.resize(p.B);
  }
};

template <typename T>
struct VertexDataPair : public Pair<VertexData<T>> {
  Pair<T> operator[](VertexPair p) const { return { this->A[p.A], this->B[p.B] }; }   // Note: this returns a copy, not a reference!
};

struct SeamInfo {
  SeamVerticesPair vertices;

  Vector2d similarityLeft;              // (a, b) treated as a 2x2 matrix [a, -b; b, a], as in Equation (2)
  Vector2d similarityRight;

  // MOSEK variables
  Variable::t similarityLeftVar;
  Variable::t similarityRightVar;
};

struct Model : public ModelBase {
  ManifoldSurfaceMesh* mesh;
  VertexPositionGeometry* geom;

  int nV;

  MatrixXd V;
  MatrixXi F;
  MatrixXd V_uv;

  MatrixXd V_uv_best;

  std::vector<util::MapDifferential> mapDiff;
  std::vector<Matrix2d> R;
  std::vector<double> c;

  // MOSEK variables
  Variable::t uv;
  Variable::t y;
};

volatile std::sig_atomic_t leaveLoop = 0;

}

// Local per-triangle computation of map differential & best fit rotation
void gmap_local(const MatrixXd& V, const MatrixXi& F, const MatrixXd& V_uv, std::vector<util::MapDifferential>& mapDiff, std::vector<Matrix2d>& R, std::vector<double>& c) {
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
void gmap_global_pre(MosekModel::t M, Model& model, std::vector<SeamInfo>& seamInfos, const std::vector<VertexPair>& landmarkRing, int iter) {
  SPDLOG_INFO("Building model for {}...", model.name);

  const int nV = model.V.rows();
  const int nF = model.F.rows();
  const int nS = (int)seamInfos.size();

  ASSERT(model.V.cols() == 3);
  ASSERT(model.F.cols() == 3);
  ASSERT(model.mapDiff.size() == nF);
  ASSERT(model.R.size() == nF);
  ASSERT(model.c.size() == nF);

  const double sqrt2 = std::sqrt(2.0);

  // Main variable (UV coordinate) per vertex
  model.uv = M->variable(model.name + "_uv", 2 * nV, Domain::unbounded());

  auto get_variable_uv_vector = [&](int vi) {
    ASSERT(vi < nV);
    return Expr::vstack(
      model.uv->index(vi),
      model.uv->index(nV + vi));
  };

  // Auxiliary variable per triangle
  Variable::t rB = M->variable(model.name + "_rB", nF, Domain::unbounded());
  Variable::t rC = M->variable(model.name + "_rC", nF, Domain::unbounded());
  Variable::t u  = M->variable(model.name + "_u" , nF, Domain::unbounded());
  Variable::t s  = M->variable(model.name + "_s" , nF, Domain::unbounded());
  Variable::t S  = M->variable(model.name + "_S" , nF, Domain::unbounded());
  Variable::t U  = M->variable(model.name + "_U" , nF, Domain::unbounded());

  std::vector<double> sqrt_area(nF);

  for (int i = 0; i < nF; ++i) {
    auto rB_i = rB->index(i);
    auto rC_i = rC->index(i);
    auto u_i  = u ->index(i);
    auto s_i  = s ->index(i);
    auto S_i  = S ->index(i);
    auto U_i  = U ->index(i);

    // Skip degenerate or inverted triangles
    if (model.mapDiff[i].Q.determinant() <= 0) {
      SPDLOG_WARN("Face {} is degenerate or inverted! Skipping", i);
      M->constraint(model.name + "_degenerate_" + std::to_string(i),
        U_i,
        Domain::equalsTo(0.0));
      continue;
    }

    auto Pinv_MOSEK = mosek_util::from_eigen_to_mosek(model.mapDiff[i].Pinv);

    sqrt_area[i] = std::sqrt(0.5 * model.mapDiff[i].P.determinant());

    // Variable 2D vectors representing the triangle corners
    int v0 = model.F(i, 0);
    int v1 = model.F(i, 1);
    int v2 = model.F(i, 2);

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
    M->constraint(model.name + "_qc_10a_B_" + std::to_string(i),
      Expr::vstack(rB_i, Expr::flatten(B)),
      Domain::inQCone());

    // |C| <= rC
    M->constraint(model.name + "_qc_10a_C_" + std::to_string(i),
      Expr::vstack(rC_i, Expr::flatten(C)),
      Domain::inQCone());

    // r = (rB + rC) / sqrt(2)
    auto r = Expr::mul(1.0 / sqrt2, Expr::add(rB_i, rC_i));

    // Equation (10b) --------------------+
    //    |C| <= tr(R^T B) - sqrt(2) u    |
    //------------------------------------+
    auto R_MOSEK = MosekMatrix::dense(mosek_util::from_eigen_to_mosek(model.R[i]));
    auto RtB = Expr::mul(R_MOSEK->transpose(), B);
    auto tr_RtB = mosek_util::trace(RtB);

    M->constraint(model.name + "_qc_10b_" + std::to_string(i),
      Expr::vstack(Expr::sub(tr_RtB, Expr::mul(sqrt2, u_i)), Expr::flatten(C)),
      Domain::inQCone());

    // Equation (10c) --------------------+
    //    sqrt((u - s)^2 + 4) <= u + s    |
    //------------------------------------+
    M->constraint(model.name + "_qc_10c_" + std::to_string(i),
      Expr::vstack(Expr::add(u_i, s_i), Expr::sub(u_i, s_i), 2.0),
      Domain::inQCone());

    // Equation (10d) ------------------------+
    //    sqrt((s / c)^2 + (c * r)^2) <= S    |
    //----------------------------------------+
    M->constraint(model.name + "_qc_10d_" + std::to_string(i),
      Expr::vstack(S_i, Expr::mul(1.0 / model.c[i], s_i), Expr::mul(model.c[i], r)),
      Domain::inQCone());

    // Equation (10e) --------------------------+
    //    sqrt(S^2 + (U - 1/4)^2) <= U + 1/4    |
    //------------------------------------------+
    M->constraint(model.name + "_qc_10e_" + std::to_string(i),
      Expr::vstack(Expr::add(U_i, 0.25), S_i, Expr::sub(U_i, 0.25)),
      Domain::inQCone());
  }

  // Final SOCP constraint ---------------------+
  //    y >= a_1 * U_1^2 + ... + a_n * U_n^2    |
  //--------------------------------------------+
  model.y = M->variable(model.name + "_y", Domain::unbounded());
  M->constraint(model.name + "_qc_y",
    Expr::vstack(0.5, model.y, Expr::mulElm(U, new_array_ptr(sqrt_area))),
    Domain::inRotatedQCone());

  // Positional constraints for two landmarks to eliminate remaining DoF (translation & rotation)
  int v0 = model.select(landmarkRing[0]).getIndex();
  int v1 = model.select(landmarkRing[landmarkRing.size() / 2]).getIndex();
  SPDLOG_INFO("Fixed vertices: {}, {}", v0, v1);

  Vector2d prev_q0 = model.V_uv.row(v0);
  Vector2d prev_q1 = model.V_uv.row(v1);

  Vector2d dir = prev_q1 - prev_q0;
  dir = { -dir[1], dir[0] };
  auto q0 = get_variable_uv_vector(v0);
  auto q1 = get_variable_uv_vector(v1);

  M->constraint(model.name + "_fix_uv_0",
    q0,
    Domain::equalsTo(0.0));

  M->constraint(model.name + "_fix_uv_1",
    Expr::dot(Expr::sub(q1, q0), mosek_util::from_eigen_to_mosek(dir)),
    Domain::equalsTo(0.0));

  if (iter % 2) {
    //-------------------------------------------------------------------------------+
    // Odd iteration: fix similarity transformations, optimize seam vertex positions |
    //-------------------------------------------------------------------------------+
    for (int i = 0; i < nS; ++i) {
      const std::string model_seam_i_str = model.name + "_seam_" + std::to_string(i);

      auto T_left = mosek_util::get_similarity_matrix(seamInfos[i].similarityLeft);
      auto T_right = mosek_util::get_similarity_matrix(seamInfos[i].similarityRight);

      SeamVertices& seamVertices = model.select(seamInfos[i].vertices);

      auto alpha0 = seamVertices[0].alphaVar = M->variable(model_seam_i_str + "_alpha_0", 2, Domain::unbounded());

      auto q0_left = get_variable_uv_vector(seamVertices[0].left.getIndex());
      auto q0_right = get_variable_uv_vector(seamVertices[0].right.getIndex());

      for (size_t j = 1; j < seamVertices.size(); ++j) {
        auto alphaj = seamVertices[j].alphaVar = M->variable(model_seam_i_str + "_alpha_" + std::to_string(j), 2, Domain::unbounded());

        auto qj_left = get_variable_uv_vector(seamVertices[j].left.getIndex());
        auto qj_right = get_variable_uv_vector(seamVertices[j].right.getIndex());

        M->constraint(model_seam_i_str + "_left_" + std::to_string(j),
          Expr::sub(Expr::sub(qj_left, q0_left), Expr::mul(T_left, Expr::sub(alphaj, alpha0))),
          Domain::equalsTo(0.0));

        M->constraint(model_seam_i_str + "_right_" + std::to_string(j),
          Expr::sub(Expr::sub(qj_right, q0_right), Expr::mul(T_right, Expr::sub(alphaj, alpha0))),
          Domain::equalsTo(0.0));
      }

      // Fix one endpoint of auxiliary seam vertex to (0, 0)
      M->constraint(model_seam_i_str + "_fix",
        alpha0,
        Domain::equalsTo(0.0));
    }

  } else {
    //------------------------------------------------------------------------------------------+
    // Even iteration: fix auxiliary seam vertex positions, optimize similarity transformations |
    //------------------------------------------------------------------------------------------+
    for (int i = 0; i < nS; ++i) {
      const std::string seam_i_str = "seam_" + std::to_string(i);
      const std::string model_seam_i_str = model.name + "_" + seam_i_str;

      // The similarity variables are shared among A & B
      if (model.name == "A") {
        seamInfos[i].similarityLeftVar = M->variable(seam_i_str + "_similarityLeft", 2, Domain::unbounded());
        seamInfos[i].similarityRightVar = M->variable(seam_i_str + "_similarityRight", 2, Domain::unbounded());
      }

      auto T_left = mosek_util::get_variable_similarity_matrix(seamInfos[i].similarityLeftVar);
      auto T_right = mosek_util::get_variable_similarity_matrix(seamInfos[i].similarityRightVar);

      SeamVertices& seamVertices = model.select(seamInfos[i].vertices);

      Vector2d alpha0 = seamVertices[0].alpha;

      auto q0_left = get_variable_uv_vector(seamVertices[0].left.getIndex());
      auto q0_right = get_variable_uv_vector(seamVertices[0].right.getIndex());

      for (size_t j = 1; j < seamVertices.size(); ++j) {
        Vector2d alphaj = seamVertices[j].alpha;

        auto alpha0j = mosek_util::from_eigen_to_mosek(Vector2d(alphaj - alpha0));

        auto qj_left = get_variable_uv_vector(seamVertices[j].left.getIndex());
        auto qj_right = get_variable_uv_vector(seamVertices[j].right.getIndex());

        M->constraint(model_seam_i_str + "_left_" + std::to_string(j),
          Expr::sub(Expr::sub(qj_left, q0_left), Expr::mul(T_left, alpha0j)),
          Domain::equalsTo(0.0));

        M->constraint(model_seam_i_str + "_right_" + std::to_string(j),
          Expr::sub(Expr::sub(qj_right, q0_right), Expr::mul(T_right, alpha0j)),
          Domain::equalsTo(0.0));
      }
    }
  }
}

void gmap_global_post(Model& model, std::vector<SeamInfo>& seamInfos, int iter) {
  SPDLOG_INFO("Retrieving solution values for {}", model.name);

  const int nV = model.V.rows();
  const int nF = model.F.rows();
  const int nS = (int)seamInfos.size();

  // Retrieve UV for each vertex
  ndarray<double, 1> uv_level = *(model.uv->level());
  for (int i = 0; i < nV; ++i)
    model.V_uv.row(i) << uv_level[i], uv_level[nV + i];

  // Retrieve auxiliary solution values
  if (iter % 2) {
    // Odd iteration: copy seam vertex positions
    for (int i = 0; i < nS; ++i) {
      SeamVertices& seamVertices = model.select(seamInfos[i].vertices);
      for (size_t j = 0; j < seamVertices.size(); ++j) {
        ndarray<double, 1> alphaj_level = *(seamVertices[j].alphaVar->level());
        seamVertices[j].alpha = kt84::vector_cast<2, Vector2d>(alphaj_level);
      }
    }

  } else if (model.name == "A") {
    // Even iteration: copy similarity transformations
    for (int i = 0; i < nS; ++i) {
      ndarray<double, 1> similarityLeft_level = *(seamInfos[i].similarityLeftVar->level());
      ndarray<double, 1> similarityRight_level = *(seamInfos[i].similarityRightVar->level());
      seamInfos[i].similarityLeft = kt84::vector_cast<2, Vector2d>(similarityLeft_level);
      seamInfos[i].similarityRight = kt84::vector_cast<2, Vector2d>(similarityRight_level);
    }
  }
}

void gmap_global(Model& modelA, Model& modelB, std::vector<SeamInfo>& seamInfos, const std::vector<VertexPair>& landmarkRing, int iter) {
  MosekModel::t M = new MosekModel("gmap");
  auto _M = finally([&]() { M->dispose(); });

  // Build SOCP
  gmap_global_pre(M, modelA, seamInfos, landmarkRing, iter);
  gmap_global_pre(M, modelB, seamInfos, landmarkRing, iter);

  // Constrain the landmarks on A & B to coincide
  auto get_variable_uv_vector = [&](const Model& model, int vi) {
    const int nV = model.mesh->nVertices();
    ASSERT(vi < nV);
    return Expr::vstack(
      model.uv->index(vi),
      model.uv->index(nV + vi));
  };

  for (size_t i = 0; i < landmarkRing.size(); ++i) {
    VertexPair v = landmarkRing[i];
    auto qA = get_variable_uv_vector(modelA, v.A.getIndex());
    auto qB = get_variable_uv_vector(modelB, v.B.getIndex());

    M->constraint("landmark_" + std::to_string(i),
      Expr::sub(qA, qB),
      Domain::equalsTo(0.0));
  }

  // Set the objective function to minimize (modelA.y + modelB.y)/4
  M->objective("obj", ObjectiveSense::Minimize, Expr::mul(0.25, Expr::add(modelA.y, modelB.y)));

  SPDLOG_INFO("Solving...");

  // Solve
  M->solve();
  SPDLOG_INFO("Primary objective value: " ANSI_BOLD ANSI_CYAN "{}" ANSI_RESET, M->primalObjValue());

  // Retrieve solution values
  gmap_global_post(modelA, seamInfos, iter);
  gmap_global_post(modelB, seamInfos, iter);
}

void gmap(
  ManifoldSurfaceMesh& A_mesh,
  ManifoldSurfaceMesh& B_mesh,
  VertexPositionGeometry& A_geom_3D,
  VertexPositionGeometry& B_geom_3D,
  const std::vector<std::set<IntPair>>& normalLandmarkGroups_i,
  const std::vector<std::set<IntPair>>& handleLandmarkGroups_i,
  int num_iter,
  std::unique_ptr<VertexPositionGeometry>& A_geom_2D,
  std::unique_ptr<VertexPositionGeometry>& B_geom_2D,
  bool debugOutput)
{
  SPDLOG_INFO(ANSI_BOLD ANSI_YELLOW "gmap BEGIN" ANSI_RESET);

  Model modelA = {"A", &A_mesh, &A_geom_3D };
  Model modelB = {"B", &B_mesh, &B_geom_3D };

  auto forEachModel = [&](std::function<void(Model& model)> f) {
    f(modelA);
    f(modelB);
  };

  forEachModel([&](Model& model) {
    // For interoperability with libigl
    model.V = util::get_vertex_position_matrix(*model.geom);
    model.F = model.mesh->getFaceVertexMatrix<int>();

    model.nV = model.V.rows();

    ASSERT_WITH_LOG(model.mesh->genus() == 0, "{}'s genus is wrong ({})", model.name, model.mesh->genus());
    ASSERT_WITH_LOG(model.mesh->nConnectedComponents() == 1, "{}'s #connectedComponents is wrong ({})", model.name, model.mesh->nConnectedComponents());
    ASSERT_WITH_LOG(model.mesh->nBoundaryLoops() == 1, "{}'s #boundaryLoops is wrong ({})", model.name, model.mesh->nBoundaryLoops());

    // Normalize area for better behavior in conformal mapping computation (which is initially equivalent to isometric mapping)
    VectorXd dblA;
    igl::doublearea(model.V, model.F, dblA);

    double total_area = 0.5 * dblA.sum();

    model.V *= 1.0 / std::sqrt(total_area);
  });

  const int nV_max = std::max<int>(modelA.nV, modelB.nV);

  std::set<VertexPair> allLandmarks;
  SetPair<Vertex> allLandmarksPerModel;

  // Convert integer-based pointer to element-based
  auto convert_landmarkGroups_i = [&](const std::vector<std::set<IntPair>>& landmarkGroups_i) {
    std::vector<std::set<VertexPair>> landmarkGroups(landmarkGroups_i.size());

    for (size_t i = 0; i < landmarkGroups_i.size(); ++i) {
      for (auto [iA, iB] : landmarkGroups_i[i]) {
        ASSERT_WITH_LOG(iA < modelA.nV, "A's landmark index {} is no less than {}", iA, modelA.nV);
        ASSERT_WITH_LOG(iB < modelB.nV, "B's landmark index {} is no less than {}", iB, modelB.nV);

        VertexPair v = {
          modelA.mesh->vertex(iA),
          modelB.mesh->vertex(iB)
        };

        landmarkGroups[i].insert(v);
        allLandmarks.insert(v);
        allLandmarksPerModel.insert(v);
      }
    }

    return landmarkGroups;
  };

  std::vector<std::set<VertexPair>> normalLandmarkGroups = convert_landmarkGroups_i(normalLandmarkGroups_i);
  std::vector<std::set<VertexPair>> handleLandmarkGroups = convert_landmarkGroups_i(handleLandmarkGroups_i);

  // The number of seams is half the number of all landmarks
  const int nS2 = (int)allLandmarks.size();
  const int nS = nS2 / 2;
  ASSERT_WITH_LOG(nS2 % 2 == 0, "Wrong number of landmarks ({})", nS2);

  // Walk along the boundary to extract boundary sides
  std::vector<VerticesPair> boundarySides(nS2);
  forEachModel([&](Model& model) {
    SPDLOG_INFO("Analyzing boundary sides on {}", model.name);

    // Starting landmark vertex
    Vertex vStart = model.select(*normalLandmarkGroups[0].begin());
    ASSERT_WITH_LOG(vStart.isBoundary(), "{}'s landmark {} is not on boundary", model.name, vStart);

    Halfedge heCurr = vStart.halfedge().twin().next();
    ASSERT(!heCurr.isInterior());
    ASSERT(heCurr.tailVertex() == vStart);

    for (size_t i = 0; i < nS2; ++i) {
      do {
        model.select(boundarySides[i]).push_back(heCurr.tailVertex());
        heCurr = heCurr.next();
      } while (model.select(allLandmarksPerModel).count(heCurr.tailVertex()) == 0);
    }

    ASSERT(heCurr.tailVertex() == vStart);
  });

  // Initialize V_uv via Tutte embedding with regular polygonal boundary
  forEachModel([&](Model& model) {
    SPDLOG_INFO("Initialization via Tutte embedding for {}", model.name);

    std::vector<std::vector<Vector2d>> boundarySidesPos(nS2);

    for (size_t i = 0; i < nS2; ++i) {
      double theta0 = -(2.0 * util::PI() * i) / nS2;
      double theta1 = -(2.0 * util::PI() * (i + 1)) / nS2;

      Vector2d pos0 = { std::cos(theta0), std::sin(theta0) };
      Vector2d pos1 = { std::cos(theta1), std::sin(theta1) };

      boundarySidesPos[i].resize(model.select(boundarySides[i]).size());

      for (size_t j = 0; j < boundarySidesPos[i].size(); ++j) {
        double a = j / (double)boundarySidesPos[i].size();
        boundarySidesPos[i][j] = (1.0 - a) * pos0 + a * pos1;
      }
    }

    // Flatten these two arrays

    VectorXi bnd(model.mesh->boundaryLoop(0).degree());
    MatrixXd bnd_uv(bnd.size(), 2);

    int index = 0;
    for (size_t i = 0; i < nS2; ++i) {
      for (size_t j = 0; j < boundarySidesPos[i].size(); ++j) {
        bnd(index) = model.select(boundarySides[i])[j].getIndex();
        bnd_uv.row(index) << boundarySidesPos[i][j].transpose();
        ++index;
      }
    }

    // Harmonic parameterization
    igl::harmonic(model.F, bnd, bnd_uv, 1, model.V_uv);

    // util::write_flattened_mesh(model.V_uv, model.F, "gmap" + model.name + "-init.obj");

    SPDLOG_INFO("Initial state conformal energy: " ANSI_BOLD ANSI_YELLOW "{}" ANSI_RESET, util::compute_conformal_energy(model.V, model.F, model.V_uv));
  });

  // Landmark group ID to mesh vertices
  VertexDataPair<int> vertexLandmarkID;
  forEachModel([&](Model& model) {
    model.select(vertexLandmarkID) = VertexData<int>(*model.mesh, -1);

    // Normal landmark group ID
    for (int i = 0; i < (int)normalLandmarkGroups.size(); ++i) {
      for (auto v : normalLandmarkGroups[i])
        model.select(vertexLandmarkID)[model.select(v)] = i;
    }

    // Handle landmark group ID (offset by nV_max)
    for (int i = 0; i < (int)handleLandmarkGroups.size(); ++i) {
      for (auto v : handleLandmarkGroups[i])
        model.select(vertexLandmarkID)[model.select(v)] = nV_max + i;
    }
  });

  forEachModel([&](Model& model) {
    SPDLOG_INFO("Landmarks on {} along the circle:", model.name);
    for (int i = 0; i < nS2; ++i) {
      VertexPair v = boundarySides[i].front();
      SPDLOG_INFO("  {}-th landmark vertex: {} (group ID: {})", i, model.select(v), model.select(vertexLandmarkID[v]));
    }
  });

  // Append front of the next boundary vertex array to the end of current seam vertex array, for convenience
  std::vector<VerticesPair> seamVertexArrays = boundarySides;
  std::vector<VertexPair> landmarkRing(nS2);
  for (int i = 0; i < nS2; ++i) {
    seamVertexArrays[i].push_back(seamVertexArrays[(i + 1) % nS2].front());
    landmarkRing[i] = seamVertexArrays[i].front();
  }

  // Find matching pair in seamVertexArrays
  std::vector<SeamInfo> seamInfos;
  seamInfos.reserve(nS);
  while (!seamVertexArrays.empty()) {
    for (auto i = seamVertexArrays.begin(); ; ++i) {
      auto [id0, id0B] = vertexLandmarkID[i->front()];
      auto [id1, id1B] = vertexLandmarkID[i->back()];
      ASSERT(id0 == id0B);
      ASSERT(id1 == id1B);
      ASSERT(id0 > -1);
      ASSERT(id1 > -1);

      auto addSeamInfo = [&](const VerticesPair& verticesLeft, const VerticesPair& verticesRight) {
        seamInfos.resize(seamInfos.size() + 1);
        SeamInfo& seamInfo = seamInfos.back();

        Pair<size_t> nV = verticesLeft.size();
        ASSERT(nV == verticesRight.size());

        seamInfo.vertices.resize(nV);

        forEachModel([&](Model& model) {
          for (size_t j = 0; j < model.select(nV); ++j) {
            model.select(seamInfo.vertices)[j].left = model.select(verticesLeft)[j];
            model.select(seamInfo.vertices)[j].right = model.select(verticesRight)[j];
            model.select(seamInfo.vertices)[j].alpha = { j / (model.select(nV) - 1.0), 0.0 };
          }
        });
      };

      if (id0 == id1 && (id0 - nV_max) % 2 == 0) {
        // Seam runs around a handle, starting from and ending at the same landmark group
        VerticesPair verticesRight = *i;
        seamVertexArrays.erase(i);

        // Look for the corresponding seam on the left
        i = seamVertexArrays.begin();
        for (; ; ++i) {
          auto [idFront, idFrontB] = vertexLandmarkID[i->front()];
          auto [idBack, idBackB] = vertexLandmarkID[i->back()];
          ASSERT(idFront == idFrontB);
          ASSERT(idBack == idBackB);

          if (idFront == id0 + 1 && idBack == id0 + 1) break;
        }
        ASSERT(i != seamVertexArrays.end());

        VerticesPair verticesLeft = *i;
        seamVertexArrays.erase(i);

        // Orient the left side the same as the right side
        verticesLeft.reverse();

        addSeamInfo(verticesLeft, verticesRight);

        // Move on to finding the next pair
        break;

      } else if (id0 < id1) {
        // Seam on the right found
        VerticesPair verticesRight = *i;
        seamVertexArrays.erase(i);

        // Look for the corresponding seam on the left
        i = seamVertexArrays.begin();
        for (; ; ++i) {
          auto [idFront, idFrontB] = vertexLandmarkID[i->front()];
          auto [idBack, idBackB] = vertexLandmarkID[i->back()];
          ASSERT(idFront == idFrontB);
          ASSERT(idBack == idBackB);

          if (idFront == id1 && idBack == id0) break;
        }
        ASSERT(i != seamVertexArrays.end());

        VerticesPair verticesLeft = *i;
        seamVertexArrays.erase(i);

        // Orient the left side the same as the right side
        verticesLeft.reverse();

        addSeamInfo(verticesLeft, verticesRight);

        // Move on to finding the next pair
        break;
      }
    }
  }
  // Report seam info
  for (int i = 0; i < nS; ++i) {
    SPDLOG_INFO("{}-th seam info:", i);
    forEachModel([&](Model& model) {
      SeamVertexPair vFront = seamInfos[i].vertices.front();
      SeamVertexPair vBack = seamInfos[i].vertices.back();
      SPDLOG_INFO("  On {}, left side: from {} to {}", model.name, model.select(vFront).left, model.select(vBack).left);
      SPDLOG_INFO("  On {}, right side: from {} to {}", model.name, model.select(vFront).right, model.select(vBack).right);
    });
  }

  // Set signal handler to allow leaving the optimization loop with CTRL+C
  SPDLOG_INFO(ANSI_BOLD "Starting G-mapping optimization!" ANSI_RESET);
  SPDLOG_INFO(ANSI_BOLD "Press Ctrl+C to leave the loop and move on with the current best result." ANSI_RESET);
  signal(SIGINT, [](int){
    SPDLOG_INFO(ANSI_BOLD "Received SIGINT; will leave the loop after the current iteration finishes." ANSI_RESET);
    leaveLoop = 1;

    // Reset signal handler to default
    std::signal(SIGINT, SIG_DFL);
  });

  // Iterative optimization
  double prev_conformal_energy;
  kt84::MinSelector<int> bestIter;
  for (int iter = 0; iter < num_iter && !leaveLoop; ++iter) {
    SPDLOG_INFO(ANSI_BOLD "Iteration {}..." ANSI_RESET, iter);

    forEachModel([&](Model& model) {
      SPDLOG_INFO("Local step for {}...", model.name);
      gmap_local(model.V, model.F, model.V_uv, model.mapDiff, model.R, model.c);

      if (iter == 0) model.c = std::vector<double>(model.F.rows(), 1.0);     // Initial values should be 1
    });

    SPDLOG_INFO("Global step...");
    gmap_global(modelA, modelB, seamInfos, landmarkRing, iter);

    // Normalize total area based on A
    double total_area = 0;
    for (int j = 0; j < modelA.F.rows(); ++j) {
      Vector2d q0 = modelA.V_uv.row(modelA.F(j, 0));
      Vector2d q1 = modelA.V_uv.row(modelA.F(j, 1));
      Vector2d q2 = modelA.V_uv.row(modelA.F(j, 2));
      Matrix2d Q;
      Q << q1 - q0, q2 - q0;
      total_area += 0.5 * Q.determinant();
    }
    modelA.V_uv *= 1.0 / std::sqrt(total_area);
    modelB.V_uv *= 1.0 / std::sqrt(total_area);

    Pair<double> conformal_energy;
    forEachModel([&](Model& model){
      int nInvertedFaces;
      model.select(conformal_energy) = util::compute_conformal_energy(model.V, model.F, model.V_uv, &nInvertedFaces);

      if (nInvertedFaces) {
        SPDLOG_INFO(ANSI_MAGENTA "{} has {} inverted faces in 2D; trying to untangle them..." ANSI_RESET, model.name, nInvertedFaces);

        if (debugOutput) {
          util::write_flattened_mesh(model.V_uv, model.F, "gmap" + model.name + "-" + std::to_string(iter) + "-tangled.obj");
        }

        UntangleParam param;
        param.theta = 0;      // Pure angle preservation
        param.debug = 0;
        // param.converge_threshold = std::numeric_limits<double>::infinity();   // Stop iteration as soon as all inverted triangles disappear
        param.maxiter = 200;
        param.bfgs_maxiter = 500;

        untangle(model.V, model.F, model.V_uv, &param);

        model.select(conformal_energy) = util::compute_conformal_energy(model.V, model.F, model.V_uv, &nInvertedFaces);

        if (nInvertedFaces) {
          SPDLOG_INFO(ANSI_BOLD ANSI_RED "  Untangling failed with {} inverted faces..." ANSI_RESET, nInvertedFaces);
        } else {
          SPDLOG_INFO(ANSI_MAGENTA "  Untangling succeeded!" ANSI_RESET);
        }
      }
    });

    prev_conformal_energy = conformal_energy.A + conformal_energy.B;
      util::compute_conformal_energy(modelB.V, modelB.F, modelB.V_uv);
    SPDLOG_INFO("Optimized conformal energy: " ANSI_BOLD ANSI_YELLOW "{}" ANSI_RESET, prev_conformal_energy);

    if (debugOutput) {
      util::write_flattened_mesh(modelA.V_uv, modelA.F, "gmapA-" + std::to_string(iter) + ".obj");
      util::write_flattened_mesh(modelB.V_uv, modelB.F, "gmapB-" + std::to_string(iter) + ".obj");
    }

    if (bestIter.update(prev_conformal_energy, iter)) {
      modelA.V_uv_best = modelA.V_uv;
      modelB.V_uv_best = modelB.V_uv;
    }
    SPDLOG_INFO("Best iteration found so far: " ANSI_BOLD ANSI_MAGENTA "{}" ANSI_RESET, bestIter.value);
  }

  SPDLOG_INFO(ANSI_BOLD ANSI_YELLOW "gmap END\a" ANSI_RESET);

  // Copy results to VertexPositionGeometry
  A_geom_2D.reset(new VertexPositionGeometry(*modelA.mesh, util::append_zero_col(modelA.V_uv_best)));
  B_geom_2D.reset(new VertexPositionGeometry(*modelB.mesh, util::append_zero_col(modelB.V_uv_best)));
}
