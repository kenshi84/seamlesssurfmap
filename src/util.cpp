#include "util.h"

#include <Eigen/Dense>
#include <igl/cat.h>
#include <igl/writeOBJ.h>

std::vector<std::string> util::tokenize(const std::string& s, char delim) {
  // https://stackoverflow.com/a/236803
  std::stringstream ss(s);
  std::string item;
  std::vector<std::string> res;
  while (std::getline(ss, item, delim)) {
    res.push_back(item);
  }
  return res;
}

MatrixXd util::get_vertex_position_matrix(const EmbeddedGeometryInterface& geom) {
  ASSERT(geom.mesh.isCompressed());

  MatrixXd V(geom.mesh.nVertices(), 3);

  for (Vertex v : geom.mesh.vertices()) {
    V.row(v.getIndex()) << kt84::vector_cast<3, Vector3d>(geom.vertexPositions[v]).transpose();
  }

  return V;
}

MatrixXd util::append_zero_col(const MatrixXd& V_uv) {
  MatrixXd Z = MatrixXd::Zero(V_uv.rows(), 1);
  MatrixXd V;
  igl::cat(2, V_uv, Z, V);
  return V;
}

void util::write_flattened_mesh(const MatrixXd& V_uv, const MatrixXi& F, const std::string& filename) {
  ASSERT(V_uv.cols() == 2);
  ASSERT(F.cols() == 3);

  SPDLOG_INFO("Writing to {} ...", filename);

  igl::writeOBJ(filename, append_zero_col(V_uv), F);
}

void util::write_point_cloud(const MatrixXd& V, const std::string& filename_base) {
  SPDLOG_INFO("Writing to {} ...", filename_base + ".ply");

  std::ofstream fout((filename_base + ".ply").c_str());
  fout
    << "ply\n"
    << "format ascii 1.0\n"
    << "element vertex " << V.rows() << "\n"
    << "property float x\n"
    << "property float y\n"
    << "property float z\n"
    << "end_header\n";
  for (int i = 0; i < V.rows(); ++i)
    fout << V.row(i) << "\n";
}

void util::write_point_cloud(const std::vector<Vector3d>& V_stdvec, const std::string& filename_base) {
  // Convert to Eigen matrix
  MatrixXd V = MatrixXd::Zero(V_stdvec.size(), 3);
  for (int i = 0; i < V.rows(); ++i)
    V.row(i) << kt84::vector_cast<3, RowVector3d>(V_stdvec[i]);

  write_point_cloud(V, filename_base);
}

void util::write_polyline(const MatrixXd& V, const std::string& filename_base) {
  SPDLOG_INFO("Writing to {} ...", filename_base + ".ply");

  // How to draw lines in .ply files using Meshlab?
  // https://stackoverflow.com/a/58648176
  /*
    ply
    format ascii 1.0
    element vertex 2
    property float x
    property float y
    property float z
    element edge 1
    property int vertex1
    property int vertex2
    end_header
    0 0 0
    0 0 1
    0 1
  */
  std::ofstream fout((filename_base + ".ply").c_str());
  fout
    << "ply\n"
    << "format ascii 1.0\n"
    << "element vertex " << V.rows() << "\n"
    << "property float x\n"
    << "property float y\n"
    << "property float z\n"
    << "element edge " << (V.rows() - 1) << "\n"
    << "property int vertex1\n"
    << "property int vertex2\n"
    << "end_header\n";
  for (int i = 0; i < V.rows(); ++i)
    fout << V.row(i) << "\n";
  for (int i = 0; i < V.rows() - 1; ++i)
    fout << i << " " << (i + 1) << "\n";
}

void util::write_polyline(const std::vector<Vector2d>& V_stdvec, const std::string& filename_base) {
  // Convert to Eigen matrix
  MatrixXd V = MatrixXd::Zero(V_stdvec.size(), 3);
  for (int i = 0; i < V.rows(); ++i)
    V.row(i) << kt84::vector_cast<2, RowVector2d>(V_stdvec[i]), 0;

  write_polyline(V, filename_base);
}

void util::write_polyline(const std::vector<Vector3d>& V_stdvec, const std::string& filename_base) {
  // Convert to Eigen matrix
  MatrixXd V = MatrixXd::Zero(V_stdvec.size(), 3);
  for (int i = 0; i < V.rows(); ++i)
    V.row(i) << kt84::vector_cast<3, RowVector3d>(V_stdvec[i]);

  write_polyline(V, filename_base);
}

void util::write_landmarks(const std::set<VertexPair>& normalLandmarks, const std::set<VertexPair>& handleLandmarks, const std::string& filename) {
  SPDLOG_INFO("Writing {} landmarks to {} ...", normalLandmarks.size() + handleLandmarks.size(), filename);
  std::ofstream fout(filename.c_str());

  std::string separator = "";

  for (const VertexPair& v : normalLandmarks) {
    fout << separator << v;
    separator = ",";
  }

  for (const VertexPair& v : handleLandmarks) {
    fout << separator << v;
    separator = ",";
  }
}

void util::polar_decomposition(const Matrix2d& A, Matrix2d& R, Matrix2d& Y) {
  // Frank Uhlig, Explicit polar decomposition and a near-characteristic polynomial: The 2×2 case, Linear Algebra and its Applications, 1981.
  // https://core.ac.uk/download/pdf/82360162.pdf

  double detA = A.determinant();
  if (detA == 0)
    SPDLOG_WARN("A's determinant is zero");

  Matrix2d B = A + std::abs(detA) * A.transpose().inverse();

  double detB = B.determinant();
  if (detB == 0)
    SPDLOG_WARN("B's determinant is zero");

  double s = 1.0 / std::sqrt(std::abs(detB));

  R = s * B;
  Y = s * (A.transpose() * A + std::abs(detA) * Matrix2d::Identity());
}

void util::compute_singular_values(const Matrix2d& A, double& sigmaLarge, double& sigmaSmall) {
  Matrix2d B = 0.5 * (A - A.transpose() + A.trace() * Matrix2d::Identity());
  Matrix2d C = A - B;

  double rB = B.norm();
  double rC = C.norm();

  const double inv_sqrt2 = 1.0 / std::sqrt(2);

  sigmaLarge = inv_sqrt2 * (rB + rC);
  sigmaSmall = inv_sqrt2 * std::abs(rB - rC);
}

double util::get_path_length(const IntrinsicGeometryInterface& geom, const std::vector<Halfedge>& path) {
  double ret = 0;
  for (auto he : path) {
    ret += geom.edgeLengths[he.edge()];
  }
  return ret;
}

void util::cut_along_path(ManifoldSurfaceMesh& mesh, EmbeddedGeometryInterface& geom, const std::vector<Halfedge>& path, std::map<Vertex, Vertex>& vertexMap_ret) {
  vertexMap_ret = {};
  for (auto he : path) {
    Halfedge he0, he1;
    std::map<Vertex, Vertex> vertexMap;
    std::tie(he0, he1, vertexMap) = mesh.separateEdge(he);
    for (auto& p : vertexMap) {
      geom.vertexPositions[p.second] = geom.vertexPositions[p.first];
    }
    vertexMap_ret.insert(vertexMap.begin(), vertexMap.end());
  }
}

void util::represent_in_local_frame(const Vector3d& p01_xyz, const Vector3d& p02_xyz, Vector2d& p01_uv, Vector2d& p02_uv) {
  // Local 2D coordinate frame
  Vector3d ez = p01_xyz.cross(p02_xyz).normalized();
  Vector3d ex = p01_xyz.normalized();
  Vector3d ey = ez.cross(ex);

  // 2D edge vectors in the local frame
  p01_uv = {ex.dot(p01_xyz), 0.};
  p02_uv = {ex.dot(p02_xyz), ey.dot(p02_xyz)};
}

util::MapDifferential util::compute_map_differential(const MatrixXd& V, const MatrixXi& F, const MatrixXd& V_uv, int faceID) {
  ASSERT(V.cols() == 3);
  ASSERT(F.cols() == 3);
  ASSERT(V_uv.cols() == 2);
  ASSERT(V_uv.rows() == V.rows());

  MapDifferential ret;

  int v0 = F(faceID, 0);
  int v1 = F(faceID, 1);
  int v2 = F(faceID, 2);

  // Original edge vectors in 3D
  Vector3d p01_xyz = (V.row(v1) - V.row(v0)).transpose();
  Vector3d p02_xyz = (V.row(v2) - V.row(v0)).transpose();

  // 2D edge vectors in the local frame
  Vector2d p01;
  Vector2d p02;
  represent_in_local_frame(p01_xyz, p02_xyz, p01, p02);

  ret.P << p01, p02;
  ret.Pinv = ret.P.inverse();

  Vector2d q01 = (V_uv.row(v1) - V_uv.row(v0)).transpose();
  Vector2d q02 = (V_uv.row(v2) - V_uv.row(v0)).transpose();

  ret.Q << q01, q02;

  ret.A = ret.Q * ret.Pinv;

  return ret;
}

double util::compute_isometric_energy(const MatrixXd& V, const MatrixXi& F, const MatrixXd& V_uv, const std::vector<double>& c) {
  ASSERT(c.size() == F.rows());

  double ret = 0;

  for (int i = 0; i < F.rows(); ++i) {
    util::MapDifferential mapDiff = compute_map_differential(V, F, V_uv, i);

    // Not well-defined for degenerate case
    if (mapDiff.A.determinant() <= 0.0)
      return std::numeric_limits<double>::infinity();

    ASSERT(c[i] > 0);
    double sigmaLarge, sigmaSmall;
    compute_singular_values(c[i] * mapDiff.A, sigmaLarge, sigmaSmall);

    double summand = squared(squared(sigmaLarge) + squared(1.0 / sigmaSmall));

    ret += summand * 0.5 * mapDiff.P.determinant();
  }
  return 0.25 * ret;
}

double util::compute_conformal_energy(const MatrixXd& V, const MatrixXi& F, const MatrixXd& V_uv, int* nInvertedFaces) {
  double ret = 0;
  if (nInvertedFaces)
    *nInvertedFaces = 0;

  for (int i = 0; i < F.rows(); ++i) {
    util::MapDifferential mapDiff = compute_map_differential(V, F, V_uv, i);

    // Not well-defined for degenerate case
    if (mapDiff.A.determinant() <= 0.0) {
      ret = std::numeric_limits<double>::infinity();
      if (nInvertedFaces)
        ++(*nInvertedFaces);
      continue;
    }

    double sigmaLarge, sigmaSmall;
    compute_singular_values(mapDiff.A, sigmaLarge, sigmaSmall);

    double summand = squared(sigmaLarge / sigmaSmall);

    ret += summand * 0.5 * mapDiff.P.determinant();
  }
  return ret;
}

double util::compute_conformal_energy(const std::array<Vector3d, 3>& p, const std::array<Vector3d, 3>& q) {
  Vector2d p01, p02;
  represent_in_local_frame(p[1] - p[0], p[2] - p[0], p01, p02);

  Vector2d q01, q02;
  represent_in_local_frame(q[1] - q[0], q[2] - q[0], q01, q02);

  Matrix2d P, Q;

  P << p01, p02;
  Q << q01, q02;

  Matrix2d A = Q * P.inverse();

  double sigmaLarge, sigmaSmall;
  compute_singular_values(A, sigmaLarge, sigmaSmall);

  return 0.5 * squared(sigmaLarge / sigmaSmall) * P.determinant();
}

void util::restore_face_vertex_ordering(const EmbeddedGeometryInterface& geomOrig, EmbeddedGeometryInterface& geomCut) {
  ManifoldSurfaceMesh& meshOrig = reinterpret_cast<ManifoldSurfaceMesh&>(geomOrig.mesh);
  ManifoldSurfaceMesh& meshCut = reinterpret_cast<ManifoldSurfaceMesh&>(geomCut.mesh);

  const int nF = meshOrig.nFaces();
  ASSERT(meshCut.nFaces() == nF);

  size_t count = 0;
  for (int i = 0; i < nF; ++i) {
    Face fOrig = meshOrig.face(i);
    Face fCut = meshCut.face(i);

    Halfedge heOrig0 = fOrig.halfedge();

    Vector3 pOrig0 = geomOrig.vertexPositions[heOrig0.vertex()];

    kt84::MinSelector<Halfedge> heCut0;
    for (Halfedge heCut : fCut.adjacentHalfedges()) {
      Vector3 pCut = geomCut.vertexPositions[heCut.vertex()];
      heCut0.update(norm(pOrig0 - pCut), heCut);
    }
    ASSERT(heCut0.score == 0);

    if (heCut0.value != fCut.halfedge()) {
      meshCut.switchFaceHalfedge(fCut, heCut0.value);
      ++count;
    }
  }
  if (count)
    SPDLOG_INFO("Restored vertex ordering for {} faces", count);
}

bool util::is_inside_triangle(const Vector2d& p0, const Vector2d& p1, const Vector2d& p2, const Vector2d& q, Vector3d& t) {
  /*
  (1 - t[1] - t[2]) * p0 + t[1] * p1 + t[2] * p2 = q
  <-->
  [p1 - p0, p2 - p0] * [t[1]] = q - p0
                       [t[2]]
  */
  Matrix2d A;
  A << p1 - p0, p2 - p0;
  ASSERT(A.determinant() != 0.0);

  Vector2d tt = A.inverse() * (q - p0);

  t[1] = tt[0];
  t[2] = tt[1];
  t[0] = 1.0 - t[1] - t[2];

  return (t.array() >= 0.0).all();
}

bool util::check_intersection(const Vector2d& p0, const Vector2d& p1, const Vector2d& q0, const Vector2d& q1, double& s, double& t) {
  /*
  (1 - s) * p0 + s * p1 = (1 - t) * q0 + t * q1
  <-->
  [p1 - p0, -q1 + q0] * [s] = -p0 + q0
                        [t]
  */
  Matrix2d A;
  A << p1 - p0, -q1 + q0;

  if (A.determinant() == 0.0) {
    s = t = std::numeric_limits<double>::quiet_NaN();
    return false;
  }

  Vector2d st = A.inverse() * (-p0 + q0);
  s = st[0];
  t = st[1];

  return s >= 0.0 && t >= 0.0 && s <= 1.0 && t <= 1.0;
}

std::set<IntPair> util::parse_landmarks(const std::string& landmarkPairs_str) {
  std::set<IntPair> ret;

  SetPair<int> checkUnique;

  for (std::string landmarkPair_str : util::tokenize(landmarkPairs_str, ',')) {
    std::vector<std::string> landmark_str = util::tokenize(landmarkPair_str, '/');
    ASSERT_WITH_LOG(landmark_str.size() == 2, "Landmarks should be specified in pairs, but was given with {} items", landmark_str.size());

    int landmark_A = std::stoi(landmark_str[0]);
    int landmark_B = std::stoi(landmark_str[1]);

    ASSERT_WITH_LOG(!checkUnique.A.count(landmark_A), "A's landmark {} is specified multiple times", landmark_A);
    ASSERT_WITH_LOG(!checkUnique.B.count(landmark_B), "B's landmark {} is specified multiple times", landmark_B);

    checkUnique.insert({ landmark_A, landmark_B });

    ret.insert({ landmark_A, landmark_B });
  }

  return ret;
}

void util::copy_mesh(ManifoldSurfaceMesh& from_mesh, VertexPositionGeometry& from_geom, std::unique_ptr<ManifoldSurfaceMesh>& to_mesh, std::unique_ptr<VertexPositionGeometry>& to_geom) {
  to_mesh = from_mesh.copy();
  to_geom = from_geom.reinterpretTo(*to_mesh);
}

void util::copy_mesh(const MeshGeometry& from, MeshGeometry& to) {
  copy_mesh(*from.mesh, *from.geom, to.mesh, to.geom);
}

bool util::compare_orig_vs_cut(ManifoldSurfaceMesh& orig_mesh, VertexPositionGeometry& orig_geom, ManifoldSurfaceMesh& cut_mesh, VertexPositionGeometry& cut_geom) {
  ASSERT(orig_mesh.isCompressed());
  ASSERT(cut_mesh.isCompressed());

  // They must have the same number of faces
  if (orig_mesh.nFaces() != cut_mesh.nFaces()) return false;

  // For each face, check if the three corners coincide exactly
  for (int i = 0; i < orig_mesh.nFaces(); ++i) {
    Halfedge heOrig = orig_mesh.face(i).halfedge();
    Halfedge heCut = cut_mesh.face(i).halfedge();

    for (int j = 0; j < 3; ++j) {
      Vector3 pOrig = orig_geom.vertexPositions[heOrig.vertex()];
      Vector3 pCut = cut_geom.vertexPositions[heCut.vertex()];

      if (pOrig != pCut) return false;

      heOrig = heOrig.next();
      heCut = heCut.next();
    }
  }
  return true;
}

SurfacePoint util::convert_sp(const SurfacePoint& spFrom, ManifoldSurfaceMesh& meshTo) {
  return SurfacePoint(meshTo.face(spFrom.face.getIndex()), spFrom.faceCoords);
}

Vertex util::convert_vertex(const Vertex& vFrom, ManifoldSurfaceMesh& meshTo) {
  SurfacePoint spTo = convert_sp(SurfacePoint(vFrom).inSomeFace(), meshTo);
  ASSERT(spTo.isOnVertex());
  return spTo.nearestVertex();
}

double util::solid_angle(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& p) {
  // L'Huilier's theorem
  // https://en.wikipedia.org/wiki/Solid_angle#Tetrahedron
  double theta_a = angle(b - p, c - p);
  double theta_b = angle(c - p, a - p);
  double theta_c = angle(a - p, b - p);
  double theta_s = 0.5 * (theta_a + theta_b + theta_c);

  double u_s = std::tan(0.5 * theta_s);
  double u_a = std::tan(0.5 * (theta_s - theta_a));
  double u_b = std::tan(0.5 * (theta_s - theta_b));
  double u_c = std::tan(0.5 * (theta_s - theta_c));

  double u = u_s * u_a * u_b * u_c;

  ASSERT(u >= 0.0);
  double tan_quater_omega = std::sqrt(u);

  double omega = 4.0 * std::atan(tan_quater_omega);

  if (dot(cross(a - p, b - p), c - p) < 0.0)
    omega *= -1.0;

  return omega;
}

int util::winding_number(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, const Vector3& p) {
  ASSERT(mesh.nBoundaryLoops() == 0);

  double windingNumberReal = 0.0;
  for (Face f : mesh.faces()) {
    Vertex v0 = f.halfedge().vertex();
    Vertex v1 = f.halfedge().next().vertex();
    Vertex v2 = f.halfedge().next().next().vertex();

    Vector3 a = geom.vertexPositions[v0];
    Vector3 b = geom.vertexPositions[v1];
    Vector3 c = geom.vertexPositions[v2];

    windingNumberReal += solid_angle(a, b, c, p);
  }
  windingNumberReal /= 4.0 * PI();

  int windingNumber = (int)std::round(windingNumberReal);
  ASSERT(std::abs(windingNumber - windingNumberReal) < 1.0e-5);

  return windingNumber;
}
