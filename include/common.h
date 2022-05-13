#pragma once

#include <Eigen/Core>
#include <fusion.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/embedded_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/surface_point.h>
#include <spdlog/spdlog.h>
#include <kt84/vector_cast.hh>
#include <kt84/MaxMinSelector.hh>

#include <vector>
#include <set>
#include <map>

using Eigen::Matrix2d;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4f;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::RowVector2d;
using Eigen::RowVector3d;
using Eigen::AlignedBox2d;
using Eigen::AlignedBox3d;

using mosek::fusion::Constraint;
using mosek::fusion::Domain;
using mosek::fusion::Expr;
using mosek::fusion::Expression;
using MosekMatrix = mosek::fusion::Matrix;
using MosekModel = mosek::fusion::Model;
using mosek::fusion::ObjectiveSense;
using mosek::fusion::Var;
using mosek::fusion::Variable;
using monty::finally;
using monty::ndarray;
using monty::new_array_ptr;

using geometrycentral::surface::Vertex;
using geometrycentral::surface::Halfedge;
using geometrycentral::surface::Corner;
using geometrycentral::surface::Edge;
using geometrycentral::surface::Face;
using geometrycentral::surface::BoundaryLoop;
using geometrycentral::surface::SurfacePoint;
using geometrycentral::surface::SurfacePointType;
using geometrycentral::surface::ManifoldSurfaceMesh;
using geometrycentral::surface::IntrinsicGeometryInterface;
using geometrycentral::surface::EmbeddedGeometryInterface;
using geometrycentral::surface::VertexPositionGeometry;
using geometrycentral::surface::readManifoldSurfaceMesh;
using geometrycentral::surface::readParameterizedManifoldSurfaceMesh;
using geometrycentral::surface::VertexData;
using geometrycentral::surface::FaceData;
using geometrycentral::surface::EdgeData;
using geometrycentral::surface::HalfedgeData;
using geometrycentral::surface::CornerData;
using geometrycentral::surface::BoundaryLoopData;
using geometrycentral::Vector2;
using geometrycentral::Vector3;

#define MAKE_FORMATTABLE(T) \
  template <> \
  struct fmt::formatter<T> : public formatter<std::string> { \
    template <typename FormatContext> \
    auto format(const T &value, FormatContext &ctx) -> decltype(ctx.out()) { \
      std::ostringstream oss; \
      oss << value; \
      return formatter<std::string>::format(oss.str(), ctx); \
    } \
  }

MAKE_FORMATTABLE(RowVector2d);
MAKE_FORMATTABLE(RowVector3d);
MAKE_FORMATTABLE(Vertex);
MAKE_FORMATTABLE(Vector3);

// Formatter for Halfedge, Edge, Face, SurfacePoint
template <>
struct fmt::formatter<Halfedge> : public formatter<std::string> {
  template <typename FormatContext>
  auto format(const Halfedge &he, FormatContext &ctx) -> decltype(ctx.out()) {
    std::ostringstream oss;
    oss << he << " (" << he.tailVertex() << " -> " << he.tipVertex() << ")";
    return formatter<std::string>::format(oss.str(), ctx);
  }
};

template <>
struct fmt::formatter<Edge> : public formatter<std::string> {
  template <typename FormatContext>
  auto format(const Edge &e, FormatContext &ctx) -> decltype(ctx.out()) {
    std::ostringstream oss;
    oss << e << " (" << e.halfedge().tailVertex() << " <-> " << e.halfedge().tipVertex() << ")";
    return formatter<std::string>::format(oss.str(), ctx);
  }
};

template <>
struct fmt::formatter<Face> : public formatter<std::string> {
  template <typename FormatContext>
  auto format(const Face &f, FormatContext &ctx) -> decltype(ctx.out()) {
    std::ostringstream oss;
    oss << f << " (" << f.halfedge().tailVertex() << ", " << f.halfedge().next().tailVertex() << ", " << f.halfedge().next().next().tailVertex() << ")";
    return formatter<std::string>::format(oss.str(), ctx);
  }
};

template <>
struct fmt::formatter<SurfacePoint> : public formatter<std::string> {
  template <typename FormatContext>
  auto format(const SurfacePoint &sp, FormatContext &ctx) -> decltype(ctx.out()) {
    std::ostringstream oss;
    oss << "[sp ";
    if (sp.type == SurfacePointType::Vertex)
      oss << sp.vertex;
    else if (sp.type == SurfacePointType::Edge)
      oss << sp.edge << " (" << sp.edge.halfedge().tailVertex() << " <-> " << sp.edge.halfedge().tipVertex() << "), tEdge = " << sp.tEdge;
    else
      oss << sp.face << " (" << sp.face.halfedge().tailVertex() << ", " << sp.face.halfedge().next().tailVertex() << ", " << sp.face.halfedge().next().next().tailVertex() << "), faceCoords = " << sp.faceCoords;
    oss << " ]";
    return formatter<std::string>::format(oss.str(), ctx);
  }
};

// https://gist.github.com/fnky/458719343aabd01cfb17a3a4f7296797
#define ANSI_BOLD    "\x1b[1m"
#define ANSI_ITALIC  "\x1b[3m"
#define ANSI_BLACK   "\x1b[30m"
#define ANSI_RED     "\x1b[31m"
#define ANSI_GREEN   "\x1b[32m"
#define ANSI_YELLOW  "\x1b[33m"
#define ANSI_BLUE    "\x1b[34m"
#define ANSI_MAGENTA "\x1b[35m"
#define ANSI_CYAN    "\x1b[36m"
#define ANSI_WHITE   "\x1b[37m"
#define ANSI_RESET   "\x1b[0m"

#define ASSERT_WITH_LOG(cond, ...) \
  do { \
    if (!(cond)) { \
      SPDLOG_ERROR(__VA_ARGS__); \
      std::abort(); \
    } \
  } while (false)

#define ASSERT(cond) ASSERT_WITH_LOG(cond, "Assertion failed")

struct MeshGeometry {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::unique_ptr<CornerData<Vector2>> uv;
};

template <typename T>
struct Pair {
  T A, B;
  bool operator==(const Pair<T>& rhs) const { return this->A == rhs.A && this->B == rhs.B; }
  bool operator<(const Pair<T>& rhs) const { return this->A == rhs.A ? this->B < rhs.B : this->A < rhs.A; }
};

struct ModelBase {
  std::string name;

  template <typename T> T& select(Pair<T>& p) const { return name == "A" ? p.A : p.B; }
  template <typename T> const T& select(const Pair<T>& p) const { return name == "A" ? p.A : p.B; }

  template <typename T> T& selectOther(Pair<T>& p) const { return name == "A" ? p.B : p.A; }
  template <typename T> const T& selectOther(const Pair<T>& p) const { return name == "A" ? p.B : p.A; }
};

using IntPair = Pair<int>;

struct VertexPair : public Pair<Vertex> {
  bool isBoundary() const {
    ASSERT(A.isBoundary() == B.isBoundary());
    return A.isBoundary();
  }
  IntPair getIndex() const { return { (int)A.getIndex(), (int)B.getIndex() }; }
};

template <typename T>
struct SetPair : public Pair<std::set<T>> {
  void insert(const Pair<T>& p) {
    this->A.insert(p.A);
    this->B.insert(p.B);
  }
};

namespace std {
inline string to_string(const IntPair& i) { return to_string(i.A) + "/" + to_string(i.B); }
inline string to_string(const VertexPair& v) { return to_string(v.getIndex()); }
}

inline std::ostream& operator<< (std::ostream& os, const VertexPair& v) { return os << std::to_string(v); }
inline std::ostream& operator<< (std::ostream& os, const IntPair& i) { return os << std::to_string(i); }

MAKE_FORMATTABLE(VertexPair);
MAKE_FORMATTABLE(IntPair);

extern const std::string VERSIONTAG;
