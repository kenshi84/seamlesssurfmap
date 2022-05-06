#include "common.h"

#include <random>

namespace util {

inline double squared(double x) { return x * x; }

inline double PI() { return std::acos(-1.0); }

std::vector<std::string> tokenize(const std::string& s, char delim);

MatrixXd get_vertex_position_matrix(const EmbeddedGeometryInterface& geom);

MatrixXd append_zero_col(const MatrixXd& V_uv);

void write_flattened_mesh(const MatrixXd& V_uv, const MatrixXi& F, const std::string& filename);

void write_point_cloud(const MatrixXd& V, const std::string& filename_base);

void write_point_cloud(const std::vector<Vector3d>& V, const std::string& filename_base);

void write_polyline(const MatrixXd& V, const std::string& filename_base);

void write_polyline(const std::vector<Vector2d>& V, const std::string& filename_base);

void write_polyline(const std::vector<Vector3d>& V, const std::string& filename_base);

void write_landmarks(const std::set<VertexPair>& normalLandmarks, const std::set<VertexPair>& handleLandmarks, const std::string& filename);

// Given nonsingular 2x2 matrix A, compute decomposition A = R * Y where
//  - R is orthonormal
//  - Y is positive definite
void polar_decomposition(const Matrix2d& A, Matrix2d& R, Matrix2d& Y);

void compute_singular_values(const Matrix2d& A, double& sigmaLarge, double& sigmaSmall);

double get_path_length(const IntrinsicGeometryInterface& geom, const std::vector<Halfedge>& path);

void cut_along_path(ManifoldSurfaceMesh& mesh, EmbeddedGeometryInterface& geom, const std::vector<Halfedge>& path, std::map<Vertex, Vertex>& vertexMap);

void represent_in_local_frame(const Vector3d& p01_xyz, const Vector3d& p02_xyz, Vector2d& p01_uv, Vector2d& p02_uv);

// Computes mapping differential A = [q01, q02] * [p01, p02]^-1 where
//   - [p01, p02] are 2D vectors corresponding to the original 3D edge vectors in V, expressed in the local coordinate frame
//   - [q01, q02] are 2D vectors corresponding to the 2D edge vectors in V_uv
struct MapDifferential {
  Matrix2d P, Pinv, Q, A;
};
MapDifferential compute_map_differential(const MatrixXd& V, const MatrixXi& F, const MatrixXd& V_uv, int faceID);

double compute_isometric_energy(const MatrixXd& V, const MatrixXi& F, const MatrixXd& V_uv, const std::vector<double>& c);

inline double compute_isometric_energy(const MatrixXd& V, const MatrixXi& F, const MatrixXd& V_uv) { return compute_isometric_energy(V, F, V_uv, std::vector<double>(F.rows(), 1.0)); }

double compute_conformal_energy(const MatrixXd& V, const MatrixXi& F, const MatrixXd& V_uv, int* nInvertedFaces = nullptr);

double compute_conformal_energy(const std::array<Vector3d, 3>& p, const std::array<Vector3d, 3>& q);

void restore_face_vertex_ordering(const EmbeddedGeometryInterface& geomOrig, EmbeddedGeometryInterface& geomCut);

// Returns true when q lies within a triangle formed by (p0, p1, p2) with barycentric coordinate t such that
// q = t[0] * p0 + t[1] * p1 + t[2] * p2,   t[0] + t[1] + t[2] = 1, t[i] >= 0
bool is_inside_triangle(const Vector2d& p0, const Vector2d& p1, const Vector2d& p2, const Vector2d& q, Vector3d& t);

// Returns true when two line segments [p0, p1] and [q0, q1] intersects at a point
// (1 - s) * p0 + s * p1 = (1 - t) * q0 + t * q1
bool check_intersection(const Vector2d& p0, const Vector2d& p1, const Vector2d& q0, const Vector2d& q1, double& s, double& t);

std::set<IntPair> parse_landmarks(const std::string& str);

void copy_mesh(ManifoldSurfaceMesh& from_mesh, VertexPositionGeometry& from_geom, std::unique_ptr<ManifoldSurfaceMesh>& to_mesh, std::unique_ptr<VertexPositionGeometry>& to_geom);

void copy_mesh(const MeshGeometry& from, MeshGeometry& to);

template <typename T>
inline Vector3 get_random_color(const T& obj, double offset = 0.3, size_t seed = 0) {
  std::hash<T> h;
  std::mt19937 gen(h(obj) + seed);
  std::uniform_real_distribution<> dist(0.0, 1.0);
  return (1.0 - offset) * Vector3{ dist(gen), dist(gen), dist(gen) } + Vector3::constant(offset);
}

// Check if cut_mesh/geom was obtained by performing cuts on orig_mesh/geom
bool compare_orig_vs_cut(ManifoldSurfaceMesh& orig_mesh, VertexPositionGeometry& orig_geom, ManifoldSurfaceMesh& cut_mesh, VertexPositionGeometry& cut_geom);
inline bool compare_orig_vs_cut(const MeshGeometry& orig, const MeshGeometry& cut) {
  return compare_orig_vs_cut(*orig.mesh, *orig.geom, *cut.mesh, *cut.geom);
}

// Assuming identical two meshes up to cutting, convert element from one mesh to the other
SurfacePoint convert_sp(const SurfacePoint& spFrom, ManifoldSurfaceMesh& meshTo);
Vertex convert_vertex(const Vertex& vFrom, ManifoldSurfaceMesh& meshTo);

template <typename T>
inline void prepend_vector(std::vector<T>& a, const std::vector<T>& b) {
  a.insert(a.begin(), b.begin(), b.end());
}

template <typename T>
inline void append_vector(std::vector<T>& a, const std::vector<T>& b) {
  a.insert(a.end(), b.begin(), b.end());
}

// Solid angle subtended by triangle (a,b,c) viewed from p
// (Not using igl::solid_angle because it is sometimes faulty somehow)
double solid_angle(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& p);

int winding_number(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, const Vector3& p);

}
