#include <geometrycentral/surface/mesh_graph_algorithms.h>

#include "lift.h"
#include "util.h"
#include "transformation.h"

namespace {

const double twoPI = 2.0 * util::PI();

// Shifts angleToShift by 2*PI so that it becomes the smallest among equivalent angles larger than angleBase (returns applied offset)
double shiftAngleAbove(double& angleToShift, double angleBase) {
  double offset = 0;
  // Shift down until it becomes smaller than base
  while (angleToShift > angleBase) {
    angleToShift -= twoPI;
    offset -= twoPI;
  }
  // Shift up untill it becomes larger than base
  while (angleToShift < angleBase) {
    angleToShift += twoPI;
    offset += twoPI;
  }
  ASSERT(angleToShift > angleBase);
  return offset;
}

// Shifts angleToShift by 2*PI so that it becomes the largest among equivalent angles smaller than angleBase (returns applied offset)
double shiftAngleBelow(double& angleToShift, double angleBase) {
  double offset = 0;
  // Shift up untill it becomes larger than base
  while (angleToShift < angleBase) {
    angleToShift += twoPI;
    offset += twoPI;
  }
  // Shift down until it becomes smaller than base
  while (angleToShift > angleBase) {
    angleToShift -= twoPI;
    offset -= twoPI;
  }
  ASSERT(angleToShift < angleBase);
  return offset;
}

bool checkAdjacent(Edge e0, Edge e1) {
  ASSERT_WITH_LOG(e0.getMesh() && e0.getMesh() == e1.getMesh(), "checkAdjacent should be called on elements from same mesh");
  ASSERT_WITH_LOG(e0 != e1, "checkAdjacent should be called on distinct edges");

  std::array<Vertex, 2> v0 = e0.adjacentVertices();
  std::array<Vertex, 2> v1 = e1.adjacentVertices();
  if (v0[0] == v1[0]) return true;
  if (v0[0] == v1[1]) return true;
  if (v0[1] == v1[0]) return true;
  if (v0[1] == v1[1]) return true;
  return false;
}

using Path = std::vector<Halfedge>;
using VertexMap = std::map<Vertex, Vertex>;
using BoundaryLoopMap = std::map<BoundaryLoop, BoundaryLoop>;

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
};

struct SeamVertexPair : public Pair<SeamVertex> {
  VertexPair left() const { return { A.left, B.left }; }
  VertexPair right() const { return { A.right, B.right }; }
};

using SeamVertices = std::vector<SeamVertex>;

struct SeamVerticesPair : public Pair<SeamVertices> {
  SeamVertexPair front() const { return { A.front(), B.front() }; }
  SeamVertexPair back() const { return { A.back(), B.back() }; }
  Pair<size_t> size() const { return { A.size(), B.size() }; }
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
};

// Data involved in analyzing G-mapping's boundary loop
struct BoundaryLoopInfo {
  size_t N;                           // Loop degree
  std::vector<Edge> edgeByIndex;      // Edges ordered counter-clockwise starting with bl.halfedge().edge()
  std::map<Edge, size_t> edgeToIndex;

  SurfacePoint param_to_edgePoint(double tLoop) const {
    while (tLoop < 0.0) tLoop += N;
    tLoop = std::fmod(tLoop, N);
    int index = (int)std::floor(tLoop);
    double tEdge = tLoop - index;
    return SurfacePoint(edgeByIndex.at(index), tEdge);
  }

  double edgePoint_to_param(Edge e, double tEdge) const {
    ASSERT(0.0 <= tEdge && tEdge < 1.0);
    return edgeToIndex.at(e) + tEdge;
  }

  std::map<double, double> intersections;     // Mapping between parameter values at intersections, always satisfying
                                              // intersections.at(intersections.at(t)) == t

  // When the two endpoints of an interval map to each other, and the winding number evaluated at a point
  // within the interval is 2, we call this interval "double-winded", which happens occasionally near
  // landmarks with high concavity.
  // For each interval (tBegin, tEnd), assume tBegin < tEnd. If the interval spans over the loop start (t=0),
  // tEnd will be larger than N.
  std::vector<std::array<double, 2>> doubleWindedIntervals;
};

struct Model : public ModelBase {
  const MeshGeometry* orig;                       // Original un-cut mesh in 3D
  const MeshGeometry* gmap;                       // Disk-topology 2D mesh resulting from G-mapping

  VertexData<Vertex> correspondingLandmark;

  EdgeData<Transformation> seamTransformation;    // Contains similarity transformation for each boundary edge
  EdgeData<Edge> correspondingEdge;               // The other edge of the corresponding edge pair on a seam

  VertexData<Halfedge> goodForPathStart;          // Contains incoming halfedge he from landmark which is valid (see below)
                                                  // if the edge is within the image of the other model in UV space; Halfedge() otherwise
  BoundaryLoopInfo bli;

  Vector2d getUV(Vertex v) const {
    return kt84::vector_cast<2, Vector2d>(gmap->geom->vertexPositions[v]);
  }

  Vector2d getUV(const SurfacePoint& sp) const {
    return kt84::vector_cast<2, Vector2d>(sp.interpolate(gmap->geom->vertexPositions));
  }

  Vector2 halfedgeVector(Halfedge he) const {
    return kt84::vector_cast<2, Vector2>(gmap->geom->halfedgeVector(he));
  }

  bool isInsideFace(const Vector2d& q, Face f, Vector3& faceCoords) const {
    std::array<Vector2d, 3> p;
    for (int i = 0; Vertex v : f.adjacentVertices())
      p[i++] = getUV(v);

    Vector3d t;
    bool ret = util::is_inside_triangle(p[0], p[1], p[2], q, t);

    faceCoords = kt84::vector_cast<3, Vector3>(t);
    return ret;
  }

  bool checkIntersection(const Vector2d& p0, const Vector2d& p1, Edge e, double& tEdge) const {
    Vector2d q0 = getUV(e.halfedge().tailVertex());
    Vector2d q1 = getUV(e.halfedge().tipVertex());

    double s, t;
    bool ret = util::check_intersection(p0, p1, q0, q1, s, t);

    tEdge = t;
    return ret;
  }

  bool checkIntersection(Edge e0, Edge e1, double& tEdge0, double& tEdge1) const {
    Vector2d p0 = getUV(e0.halfedge().tailVertex());
    Vector2d p1 = getUV(e0.halfedge().tipVertex());

    Vector2d q0 = getUV(e1.halfedge().tailVertex());
    Vector2d q1 = getUV(e1.halfedge().tipVertex());

    return util::check_intersection(p0, p1, q0, q1, tEdge0, tEdge1);
  }

  // Find shortest path from vSource to a landmark vertex, such that the first halfedge is within the image of the other model in UV space
  std::vector<Halfedge> pathToNearestLandmark(Vertex vSource) const {
    // If vSource is already adjacent to landmark and the halfedge is suitable, return it
    if (goodForPathStart[vSource] != Halfedge()) {
      return { goodForPathStart[vSource].twin() };
    }

    VertexData<double> dist(*gmap->mesh, std::numeric_limits<double>::infinity());
    VertexData<Halfedge> prev(*gmap->mesh);

    std::set<Vertex> candidate;
    std::set<Vertex> finalized = { vSource };

    // Initialize { dist, prev, candidate }
    dist[vSource] = 0;
    for (Halfedge heOut : vSource.outgoingHalfedges()) {
      Vertex vNeighbor = heOut.tipVertex();

      // Path should reach landmark only via goodForPathStart
      if (correspondingLandmark[vNeighbor] != Vertex()) continue;

      dist[vNeighbor] = gmap->geom->edgeLengths[heOut.edge()];
      prev[vNeighbor] = heOut;
      candidate.insert(vNeighbor);
    }

    while (true) {
      // Pick the current vertex in candidate with smallest distance
      kt84::MinSelector<Vertex> vCurr;
      for (Vertex v : candidate) {
        vCurr.update(dist[v], v);
      }

      // Move it from candidate to finalized
      candidate.erase(vCurr.value);
      finalized.insert(vCurr.value);

      // Termination condition met; trace path and return
      if (goodForPathStart[vCurr.value] != Halfedge()) {
        std::vector<Halfedge> ret = { goodForPathStart[vCurr.value].twin() };
        for (Vertex v = vCurr.value; v != vSource; ) {
          ret.push_back(prev[v]);
          v = prev[v].tailVertex();
        }
        // Make it from vSource to the nearest landmark
        std::reverse(ret.begin(), ret.end());
        return ret;
      }

      for (Halfedge heOut : vCurr.value.outgoingHalfedges()) {
        Vertex vNeighbor = heOut.tipVertex();

        // Skip if already finalized
        if (finalized.count(vNeighbor)) continue;

        // Path should reach landmark only via goodForPathStart
        if (correspondingLandmark[vNeighbor] != Vertex()) continue;

        // Dijkstra update rule
        double distNew = dist[vCurr.value] + gmap->geom->edgeLengths[heOut.edge()];
        if (distNew < dist[vNeighbor]) {
          dist[vNeighbor] = distNew;
          prev[vNeighbor] = heOut;
          candidate.insert(vNeighbor);
        }
      }

      if (candidate.empty()) {
        // Shouldn't reach here...
        return {};
      }
    }
  }

  // For a point on a boundary edge, check if it is in a double-winded interval.
  // If so, project the point onto some other face and return it as well.
  bool isBoundaryEdgePointInDoubleWindedInterval(Edge e, double tEdge, SurfacePoint& facePoint) const {
    ASSERT(e.isBoundary());

    double tLoop = bli.edgePoint_to_param(e, tEdge);

    facePoint = {};

    for (auto [tBegin, tEnd] : bli.doubleWindedIntervals) {
      // Check if tLoop is within the interval (check also the one shifted by N)
      if ((tBegin < tLoop && tLoop < tEnd) || (tBegin < (tLoop + bli.N) && (tLoop + bli.N) < tEnd)) {
        // Evaluate UV of the edge point
        Vector2d p = getUV(SurfacePoint(e, tEdge));

        // Find face containing p (and not adjacent to e)
        for (Face f : gmap->mesh->faces()) {
          if (f != e.halfedge().face()) {
            Vector3 faceCoords;
            if (isInsideFace(p, f, faceCoords)) {
              // There should be just one such face
              if (facePoint != SurfacePoint()) {
                SPDLOG_WARN("  Strange situation detected: edge point {} corresponds to more than one face points, {} and {}",
                  SurfacePoint(e, tEdge), facePoint, SurfacePoint(f, faceCoords));
              }

              facePoint = SurfacePoint(f, faceCoords);
            }
          }
        }

        // One face should have been found
        ASSERT(facePoint != SurfacePoint());
        return true;
      }
    }
    return false;
  }

  // FOR DEBUGGING ONLY
  std::vector<SurfacePoint> findCorrespondingPointFromUV(const Vector2d& uv) const {
    std::vector<SurfacePoint> res;
    for (Face f : gmap->mesh->faces()) {
      Vector3 faceCoords;
      if (isInsideFace(uv, f, faceCoords)) {
        res.emplace_back(f, faceCoords);
      }
    }
    return res;
  }

  // FOR DEBUGGING ONLY
  std::vector<Vector3d> easyLiftPathPoints(const std::vector<Vector2d>& pathPoints) const {
    std::vector<Vector3d> res;
    for (const Vector2d& uv : pathPoints) {
      auto sp_cut = findCorrespondingPointFromUV(uv);
      if (sp_cut.size() == 1) {
        SurfacePoint sp_orig = util::convert_sp(sp_cut[0], *orig->mesh);
        res.push_back(kt84::vector_cast<3, Vector3d>(sp_orig.interpolate(orig->geom->vertexPositions)));
      }
    }
    return res;
  }
};

// Declare models as global
Model modelA;
Model modelB;

void forEachModel(std::function<void(Model& model)> f) {
  f(modelA);
  f(modelB);
}

void forEachModel2(std::function<void(Model& modelP, Model& modelQ)> f) {
  f(modelA, modelB);
  f(modelB, modelA);
}

} // anonymous namespace end

int count = 0;

std::vector<std::vector<SurfacePoint>> lift_main_recurse(const Model& modelQ, Vertex vQ0, std::vector<SurfacePoint> oneRing, std::vector<Vector2d> pathPoints) {
  std::vector<SurfacePoint> pathQ = { oneRing[0] };       // pathQ includes the starting point

  const size_t N = pathPoints.size() - 1;

  for (size_t i = 0; i < N; ++i) {
    // Every point in pathPoints except pathPoints[0] must be found in some face of Q
    bool foundInFace = false;

    // Repeat until such a face is found
    while (true) {
      // Consider the first segment
      const Vector2d& p0 = pathPoints[0];
      const Vector2d& p1 = pathPoints[1];

      // Check if p1 is within one-ring
      for (const SurfacePoint& fpQ /* face point in Q */ : oneRing) {
        Vector3 faceCoords;
        if (modelQ.isInsideFace(p1, fpQ.face, faceCoords)) {
          foundInFace = true;

          // Add it to pathQ
          pathQ.push_back(SurfacePoint(fpQ.face, faceCoords));

          // Update oneRing
          oneRing = { pathQ.back() };

          // Pop front of pathPoints
          pathPoints.erase(pathPoints.begin());

          break;
        }
      }
      if (foundInFace) break;

      // p1 is not within one-ring
      // Look for an intersection between [p0, p1] and the rim of one-ring
      bool foundOnEdge = false;
      for (const SurfacePoint& fpQ : oneRing) {
        // The number of edges to test depends on the type of one-ring
        std::set<Edge> edgesToTest;

        if (fpQ.isOnVertex()) {
          // Special case: fpQ represents vQ0 --> test against the opposite edge
          ASSERT(fpQ.nearestVertex() == vQ0);

          for (Halfedge heQ : fpQ.face.adjacentHalfedges()) {
            if (heQ.tailVertex() != vQ0 && heQ.tipVertex() != vQ0) {
              edgesToTest.insert(heQ.edge());
            }
          }
          ASSERT(edgesToTest.size() == 1);

        } else if (fpQ.isOnEdge()) {
          // The common case: fpQ is on an edge of fpQ.face --> test against the other two edges
          for (Edge eQ : fpQ.face.adjacentEdges()) {
            if (eQ != fpQ.nearestEdge())
              edgesToTest.insert(eQ);
          }
          ASSERT(edgesToTest.size() == 2);

        } else {
          // The less common case: fpQ is inside fpQ.face, corresponding to a vertex in pathP --> test against all three edges of fpQ.face
          for (Edge eQ : fpQ.face.adjacentEdges()) {
            edgesToTest.insert(eQ);
          }
          ASSERT(edgesToTest.size() == 3);
        }

        // Check for intersection between the segment [p0, p1] and each item in edgesToTest
        for (Edge eQ : edgesToTest) {
          double tEdge;
          if (modelQ.checkIntersection(p0, p1, eQ, tEdge)) {
            foundOnEdge = true;

            // Update pathQ and pathPoints[0], regardless of situation
            pathQ.push_back(SurfacePoint(eQ, tEdge).inSomeFace());
            pathPoints[0] = modelQ.getUV(pathQ.back());

            bool isBranching = false;

            std::vector<std::vector<SurfacePoint>> branchedPathQ;

            // If the edge is on boundary and the edge point belongs to a double-winded interval, project the point
            // onto some face containing it (one must exist), and branch to there
            SurfacePoint fpQ_projected;
            if (eQ.isBoundary() && modelQ.isBoundaryEdgePointInDoubleWindedInterval(eQ, tEdge, fpQ_projected)) {
              isBranching = true;

              // Call recursively
              util::append_vector(branchedPathQ, lift_main_recurse(modelQ, vQ0, { fpQ_projected }, pathPoints));
            }

            if (eQ.isBoundary()) {
              // Apply seam transformation to pathPoints
              for (Vector2d& p : pathPoints) {
                p = modelQ.seamTransformation[eQ].apply(p);
              }

#if 0
              util::write_polyline(pathPoints, "pathPoints-" + std::to_string(count));
              util::write_point_cloud(modelQ.easyLiftPathPoints(pathPoints), "pathPoints-" + std::to_string(count) + "-" + modelQ.name + "-orig");
              ++count;
#endif

              Edge eQNext = modelQ.correspondingEdge[eQ];
              double tEdgeNext = 1.0 - tEdge;               // The other edge is oriented the other way

              // If the seam-transformed boundary edge point belongs to a double-winded interval, project the
              // point onto some face containing it (one must exist), and branch to there
              if (modelQ.isBoundaryEdgePointInDoubleWindedInterval(eQNext, tEdgeNext, fpQ_projected)) {
                isBranching = true;

                // Call recursively
                util::append_vector(branchedPathQ, lift_main_recurse(modelQ, vQ0, { fpQ_projected }, pathPoints));
              }

              if (isBranching) {
                // Call recursively
                util::append_vector(branchedPathQ, lift_main_recurse(modelQ, vQ0, { SurfacePoint(eQNext, tEdgeNext).inSomeFace() }, pathPoints));

                // Prepend upstream pathQ to each branchedPathQ (after erasing the front of each branchedPathQ)
                for (size_t k = 0; k < branchedPathQ.size(); ++k) {
                  branchedPathQ[k].erase(branchedPathQ[k].begin());
                  util::prepend_vector(branchedPathQ[k], pathQ);
                }

                return branchedPathQ;
              }

              // Update oneRing
              oneRing = { SurfacePoint(eQNext, tEdgeNext).inSomeFace() };

            } else {
              // Simple base case of crossing an interior edge;

              // Get face on the other side of eQ
              Face fQNext = eQ.halfedge().face() == fpQ.face ? eQ.halfedge().twin().face() : eQ.halfedge().face();

              // Update oneRing as edge point of that face
              oneRing = { SurfacePoint(eQ, tEdge).inFace(fQNext) };
            }

            break;
          }
        }

        if (foundOnEdge) break;
      }
      // There must exist an intersection
      ASSERT(foundOnEdge);
    } // Continue to the next segment, until foundInFace becomes true
    // The loop finishes when p1 is found in a face
    ASSERT(foundInFace);
  }

  return { pathQ };
}

void lift_main_root(const Model& modelP, const Model& modelQ, Vertex vP_target, std::vector<Halfedge>& pathP, std::vector<std::vector<SurfacePoint>>& pathQ) {
  // If vP_target is a landmark, return its corresponding landmark
  if (modelP.correspondingLandmark[vP_target] != Vertex()) {
    pathQ = {{ SurfacePoint(modelP.correspondingLandmark[vP_target]).inSomeFace() }};
    return;
  }

  // Find path to the nearest landmark
  pathP = modelP.pathToNearestLandmark(vP_target);

  // Reverse the path direction (from landmark to x)
  std::reverse(pathP.begin(), pathP.end());
  for (size_t i = 0; i < pathP.size(); ++i) {
    pathP[i] = pathP[i].twin();
  }

  Vertex vP0 = pathP[0].tailVertex();
  Vertex vQ0 = modelP.correspondingLandmark[vP0];

  // Initialize one-ring on Q
  std::vector<SurfacePoint> oneRing;                      // The size of one-ring is initially the number of vQ0's adjacent faces
  for (Face fQ : vQ0.adjacentFaces()) {                   // While lifting, the size is always 1 (the face on which the found point lies)
    oneRing.push_back(SurfacePoint(vQ0).inFace(fQ));
  }

  // Initialize path points in UV space
  std::vector<Vector2d> pathPoints;
  pathPoints.reserve(pathP.size() + 1);
  for (size_t i = 0; i < pathP.size(); ++i) {
    if (i == 0) {
      pathPoints.emplace_back(modelP.getUV(pathP[i].tailVertex()));
    }
    pathPoints.emplace_back(modelP.getUV(pathP[i].tipVertex()));
  }

#if 0
  util::write_polyline(pathPoints, "pathPoints-" + std::to_string(count));
  util::write_point_cloud(modelQ.easyLiftPathPoints(pathPoints), "pathPoints-" + std::to_string(count) + "-" + modelQ.name + "-orig");
  ++count;
#endif

  pathQ = lift_main_recurse(modelQ, vQ0, oneRing, pathPoints);
}

void lift_prepare(
  const MeshGeometry& A_orig,
  const MeshGeometry& B_orig,
  const MeshGeometry& A_gmap,
  const MeshGeometry& B_gmap,
  const std::vector<std::set<IntPair>>& normalLandmarkGroups_i,
  const std::vector<std::set<IntPair>>& handleLandmarkGroups_i)
{
  modelA = { "A", &A_orig, &A_gmap };
  modelB = { "B", &B_orig, &B_gmap };

  forEachModel([&](Model& model) {
    ASSERT_WITH_LOG(model.orig->mesh->nConnectedComponents() == 1, "{}-orig has multiple connected components ({})", model.name, model.orig->mesh->nConnectedComponents());
    ASSERT_WITH_LOG(model.orig->mesh->nBoundaryLoops() == 0, "{}-orig has boundary loops ({})", model.name, model.orig->mesh->nBoundaryLoops());

    ASSERT_WITH_LOG(model.gmap->mesh->genus() == 0, "{}-gmap has non-zero genus ({})", model.name, model.gmap->mesh->genus());
    ASSERT_WITH_LOG(model.gmap->mesh->nConnectedComponents() == 1, "{}-gmap has multiple connected components ({})", model.name, model.gmap->mesh->nConnectedComponents());
    ASSERT_WITH_LOG(model.gmap->mesh->nBoundaryLoops() == 1, "{}-gmap has multiple boundary loops ({})", model.name, model.gmap->mesh->nBoundaryLoops());

    model.gmap->geom->requireVertexAngleSums();
    model.gmap->geom->requireEdgeLengths();       // Needed despite having already required vertexAngleSums, because interanlly
                                                // EmbeddedGeometryInterface::computeCornerAngles() is called which computes
                                                // cornerAngles from vertexPositions, not from edgeLengths

    // Fill helper arrays with null
    model.correspondingLandmark = VertexData<Vertex>(*model.gmap->mesh, Vertex());
    model.seamTransformation = EdgeData<Transformation>(*model.gmap->mesh, Transformation());
    model.correspondingEdge = EdgeData<Edge>(*model.gmap->mesh, Edge());
    model.goodForPathStart = VertexData<Halfedge>(*model.gmap->mesh, Halfedge());
  });
  ASSERT_WITH_LOG(modelA.orig->mesh->genus() == modelB.orig->mesh->genus(), "A-orig and B-orig have different genus ({} vs {})", modelA.orig->mesh->genus(), modelB.orig->mesh->genus());

  const int nV_max = std::max<int>(
    modelA.gmap->mesh->nVertices(),
    modelB.gmap->mesh->nVertices());

  std::set<VertexPair> allLandmarks;
  SetPair<Vertex> allLandmarksPerModel;

  // Convert integer-based pointer to element-based
  auto convert_landmarkGroups_i = [&](const std::vector<std::set<IntPair>>& landmarkGroups_i) {
    std::vector<std::set<VertexPair>> landmarkGroups(landmarkGroups_i.size());

    for (size_t i = 0; i < landmarkGroups_i.size(); ++i) {
      for (auto [iA, iB] : landmarkGroups_i[i]) {
        ASSERT_WITH_LOG(iA < modelA.gmap->mesh->nVertices(), "A's landmark index {} is no less than {}", iA, modelA.gmap->mesh->nVertices());
        ASSERT_WITH_LOG(iB < modelB.gmap->mesh->nVertices(), "B's landmark index {} is no less than {}", iB, modelB.gmap->mesh->nVertices());

        VertexPair v = {
          modelA.gmap->mesh->vertex(iA),
          modelB.gmap->mesh->vertex(iB)
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
  ASSERT(nS2 % 2 == 0);

  // Walk along the boundary to extract boundary sides
  std::vector<VerticesPair> boundarySides(nS2);
  forEachModel([&](Model& model) {
    SPDLOG_INFO("Analyzing boundary sides on {}", model.name);

    // Starting landmark vertex
    Vertex vStart = model.select(*normalLandmarkGroups[0].begin());
    ASSERT(vStart.isBoundary());

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

    for (Vertex v : model.gmap->mesh->vertices()) {
      double d = model.gmap->geom->vertexAngleSums[v] - twoPI;
      if (d > 0.1) {
        SPDLOG_INFO("  Vertex {} angle sum is larger than 2*pi by {}, isLandmark = {}, isBoundary = {}", v, d, model.select(allLandmarksPerModel).count(v) > 0, v.isBoundary());
      }
    }
  });

  // Landmark group ID to mesh vertices
  VertexDataPair<int> vertexLandmarkID;
  forEachModel([&](Model& model) {
    model.select(vertexLandmarkID) = VertexData<int>(*model.gmap->mesh, -1);

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

  // Fill correspondingLandmark
  for (auto [vA, vB] : allLandmarks) {
    modelA.correspondingLandmark[vA] = vB;
    modelB.correspondingLandmark[vB] = vA;
  }

  // Fill seamTransformation & correspondingEdge
  for (const SeamInfo& seamInfo : seamInfos) {
    SeamVertexPair vFront = seamInfo.vertices.front();
    SeamVertexPair vBack = seamInfo.vertices.back();

    VertexPair vFrontLeft = vFront.left();
    VertexPair vFrontRight = vFront.right();
    VertexPair vBackLeft = vBack.left();
    VertexPair vBackRight = vBack.right();

    Vector2d pFrontLeft = modelA.getUV(vFrontLeft.A);
    Vector2d pFrontRight = modelA.getUV(vFrontRight.A);
    Vector2d pBackLeft = modelA.getUV(vBackLeft.A);
    Vector2d pBackRight = modelA.getUV(vBackRight.A);

    // Similarity transformation bringing the left side to the right side
    Transformation tfmLeft = Transformation::make(pFrontLeft, pBackLeft, pFrontRight, pBackRight);

    // The opposite one
    Transformation tfmRight = tfmLeft.inverse();

    // Assign transformation to every edge on the seam
    forEachModel([&](Model& model) {
      Halfedge heLeft = model.select(vFrontLeft).halfedge().twin();
      Halfedge heRight = model.select(vFrontRight).halfedge().twin().next();
      for (size_t i = 0; i < model.select(seamInfo.vertices.size()) - 1; ++i) {
        model.seamTransformation[heLeft.edge()] = tfmLeft;
        model.seamTransformation[heRight.edge()] = tfmRight;

        model.correspondingEdge[heLeft.edge()] = heRight.edge();
        model.correspondingEdge[heRight.edge()] = heLeft.edge();

        heLeft = heLeft.prevOrbitVertex();
        heRight = heRight.next();
      }
      ASSERT(model.select(allLandmarksPerModel).count(heLeft.tipVertex()));
      ASSERT(model.select(allLandmarksPerModel).count(heRight.tailVertex()));
    });

  }

  // Fill goodForPathStart
  Pair<bool> found;
  for (VertexPair landmark : allLandmarks) {
    // Landmark vertex can be locally non-injective, which cannot be used for lifting
    if (modelA.gmap->geom->vertexAngleSums[landmark.A] >= twoPI) continue;
    if (modelB.gmap->geom->vertexAngleSums[landmark.B] >= twoPI) continue;

    // Test I: do boundary edges adjacent to the landmark self-intersect with other boundary edges?
    Pair<bool> selfIntersected = { false, false };
    forEachModel([&](Model& model){
      // Take two boundary edges adjacent to the landmark
      Halfedge he0 = model.select(landmark).halfedge().twin();
      Halfedge he1 = he0.next();
      Edge e[2] = { he0.edge(), he1.edge() };

      // For each of these two edges, check intersection with other boundary edges
      for (int i = 0; i < 2; ++i) {
        // Get endpoints UV of e[i]
        std::array<Vertex, 2> v = e[i].adjacentVertices();
        Vector2d p0 = model.getUV(v[0]);
        Vector2d p1 = model.getUV(v[1]);

        // Test intersection with other non-adjacent boundary edges
        for (Edge eOther : model.gmap->mesh->edges()) {
          if (eOther != e[i] && !checkAdjacent(eOther, e[i]) && eOther.isBoundary()) {
            double tEdge;
            if (model.checkIntersection(p0, p1, eOther, tEdge)) {
              model.select(selfIntersected) = true;
              return;
            }
          }
        }
      }
    });
    if (selfIntersected.A || selfIntersected.B) {
      SPDLOG_INFO("Landmark {} has self-intersection nearby, so skipping", landmark);
      continue;
    }

    // Test II: check if the triangle fans of A & B intersect in a simple way

    // Collect necessary info (halgedge & angle)
    Pair<Halfedge> heBegin;
    Pair<Halfedge> heEnd;

    Pair<double> angleBegin;
    Pair<double> angleEnd;

    forEachModel([&](Model& model){
      model.select(heBegin) = model.select(landmark).halfedge();
      model.select(heEnd) = model.select(landmark).halfedge().twin().next();

      model.select(angleBegin) = arg(model.halfedgeVector(model.select(heBegin)));
      model.select(angleEnd) = arg(model.halfedgeVector(model.select(heEnd)));

      // Ensure angleBegin < angleEnd
      shiftAngleAbove(model.select(angleEnd), model.select(angleBegin));
    });

    // Register suitable heBegin and/or heEnd
    forEachModel2([&](Model& modelP, Model& modelQ){
      double P_angleBegin = modelP.select(angleBegin);
      double P_angleEnd = modelP.select(angleEnd);

      double Q_angleBegin = modelQ.select(angleBegin);
      double Q_angleEnd = modelQ.select(angleEnd);

      Halfedge P_heBegin = modelP.select(heBegin);
      Halfedge P_heEnd = modelP.select(heEnd);

      // Ensure Q_angleBegin < P_angleBegin
      P_angleEnd += shiftAngleAbove(P_angleBegin, Q_angleBegin);

      // Case where P_angleBegin is in [Q_angleBegin, Q_angleEnd]
      if (P_angleBegin < Q_angleEnd) {
        if (P_angleEnd < Q_angleEnd) {
          // The simplest case: P's fan is entirely cotained in Q's fan -> register both halfedges as good
          modelP.goodForPathStart[P_heBegin.tipVertex()] = P_heBegin;
          modelP.goodForPathStart[P_heEnd.tipVertex()] = P_heEnd;
          modelP.select(found) = true;

          SPDLOG_INFO("{}'s goodForPathStart: {}", modelP.name, P_heBegin.tipVertex());
          SPDLOG_INFO("{}'s goodForPathStart: {}", modelP.name, P_heEnd.tipVertex());

        } else if (P_angleEnd < Q_angleBegin + twoPI) {
          // Another common case: the two fans overlap in a simple way, producing just one intersected fan -> register P_heBegin as good
          modelP.goodForPathStart[P_heBegin.tipVertex()] = P_heBegin;
          modelP.select(found) = true;

          SPDLOG_INFO("{}'s goodForPathStart: {}", modelP.name, P_heBegin.tipVertex());

        } else {
          // Relatively rare case: the two fans overlap in a differnt way, producing two intersected fans -> ignore
        }

      // Otherwise, make P_angleBegin below Q_angleBegin; if the two fans overlap in a simple way, register P_heEnd as good
      } else {
        // Ensure P_angleBegin < Q_angleBegin
        P_angleEnd += shiftAngleBelow(P_angleBegin, Q_angleBegin);

        if (Q_angleBegin < P_angleEnd && P_angleEnd < Q_angleEnd) {
          modelP.goodForPathStart[P_heEnd.tipVertex()] = P_heEnd;
          modelP.select(found) = true;

          SPDLOG_INFO("{}'s goodForPathStart: {}", modelP.name, P_heEnd.tipVertex());
        }
      }
    });
  }

  ASSERT_WITH_LOG(found.A, "goodForPathStart is empty on A");
  ASSERT_WITH_LOG(found.B, "goodForPathStart is empty on B");

  // Fill BoundaryLoopInfo
  forEachModel([&](Model& model){
    SPDLOG_INFO("Analyzing boundary loop on {}...", model.name);

    BoundaryLoop bl = model.gmap->mesh->boundaryLoop(0);

    BoundaryLoopInfo& bli = model.bli;

    // Cache degree
    bli.N = bl.degree();

    // Fill edgeByIndex
    bli.edgeByIndex.reserve(bli.N);
    for (Halfedge he : bl.adjacentHalfedges()) {    // This goes clockwise, so we must reverse the order afterwards
      bli.edgeByIndex.push_back(he.edge());
    }
    std::reverse(bli.edgeByIndex.begin(), bli.edgeByIndex.end());

    // Fill edgeToIndex
    for (size_t i = 0; i < bli.N; ++i) {
      bli.edgeToIndex[bli.edgeByIndex.at(i)] = i;
    }

    // Compute intersections among non-adjacent edge pairs
    for (size_t i0 = 0; i0 < bli.N; ++i0) {
      for (size_t i1 = i0 + 1; i1 < bli.N; ++i1) {
        // Skip adjacent pairs
        if (i1 == (i0 + 1) % bli.N) continue;
        if (i1 == (i0 + bli.N - 1) % bli.N) continue;

        Edge e0 = bli.edgeByIndex.at(i0);
        Edge e1 = bli.edgeByIndex.at(i1);
        double tEdge0;
        double tEdge1;
        if (model.checkIntersection(e0, e1, tEdge0, tEdge1)) {
          double tLoop0 = bli.edgePoint_to_param(e0, tEdge0);
          double tLoop1 = bli.edgePoint_to_param(e1, tEdge1);

          bli.intersections[tLoop0] = tLoop1;
          bli.intersections[tLoop1] = tLoop0;
        }
      }
    }

    // Fill doubleWindedIntervals
    for (auto i0 = bli.intersections.begin(); i0 != bli.intersections.end(); ++i0) {
      auto i1 = i0;

      bool spansLoopStart = i0 == --bli.intersections.end();
      if (spansLoopStart) {
        i1 = bli.intersections.begin();
      } else {
        ++i1;
      }

      // Check if these two intersections map to each other
      if (i0->second == i1->first) {
        ASSERT(i1->second == i0->first);

        double t0 = i0->first;
        double t1 = i1->first;

        if (spansLoopStart)
          t1 += bli.N;

        // Evaluate UV of midpoint of the interval
        double tMid = 0.5 * (t0 + t1);

        SurfacePoint edgePoint = bli.param_to_edgePoint(tMid);

        Vector2 p = kt84::vector_cast<2, Vector2>(model.getUV(edgePoint));

        // Compute winding number at p by summing up signed angles subtended by other boundary edges
        double windingNumberReal = 0.0;
        for (Edge eOther : bli.edgeByIndex) {
          if (eOther != edgePoint.edge) {
            Vector2 q = kt84::vector_cast<2, Vector2>(model.getUV(eOther.halfedge().tailVertex()));
            Vector2 r = kt84::vector_cast<2, Vector2>(model.getUV(eOther.halfedge().tipVertex()));

            windingNumberReal += orientedAngle(q - p, r - p);
          }
        }

        // Add contribution from edgePoint.edge
        windingNumberReal += util::PI();

        // Normalize to integral value
        windingNumberReal /= twoPI;

        // Round to integer
        int windingNumber = (int)std::round(windingNumberReal);
        ASSERT(std::abs(windingNumberReal - windingNumber) < 1e-5);

        ASSERT(windingNumber >= 1);

        if (windingNumber > 2)
          SPDLOG_WARN("  Edge point ({}, tEdge = {}) has winding number {}, very weird", edgePoint.edge, edgePoint.tEdge, windingNumber);

        // Double-winded interval found
        if (windingNumber == 2) {
#if 0
          // Write it as polyline
          std::vector<Vector2d> points;
          for (double tLoop = t0; tLoop <= t1; ) {
            points.push_back(model.getUV(bli.param_to_edgePoint(tLoop)));

            // Increase tLoop depending on its current value
            if (tLoop == t0 || tLoop == t1) {
              tLoop = std::ceil(tLoop);
            } else if (tLoop + 1 > t1) {
              tLoop = t1;
            } else {
              ++tLoop;
            }
          }
          util::write_polyline(points, "doubleWindedInterval-" + model.name + "-" + std::to_string(bli.doubleWindedIntervals.size()));
#endif

          SPDLOG_INFO("  Boundary loop interval containing ({}, tEdge = {}) is double-winded", edgePoint.edge, edgePoint.tEdge);
          bli.doubleWindedIntervals.push_back({t0, t1});
        }
      }
    }

    // Sanity check that the number of double-winded intervals is less than or equal to the turning number minus 1
    double turningNumberReal = 0.0;
    for (size_t i = 0; i < bli.N; ++i) {
      Edge e0 = bli.edgeByIndex.at(i);
      Edge e1 = bli.edgeByIndex.at(i == bli.N - 1 ? 0 : i + 1);

      Vector2 u = model.halfedgeVector(e0.halfedge());
      Vector2 v = model.halfedgeVector(e1.halfedge());

      turningNumberReal += orientedAngle(u, v);
    }
    turningNumberReal /= twoPI;

    int turningNumber = (int)std::round(turningNumberReal);
    ASSERT(std::abs(turningNumberReal - turningNumber) < 1e-5);

    ASSERT(turningNumber >= 1);

    ASSERT((int)bli.doubleWindedIntervals.size() <= turningNumber - 1);
  });
}


void lift(
  const MeshGeometry& A_orig, const MeshGeometry& B_orig,
  const MeshGeometry& A_gmap, const MeshGeometry& B_gmap,
  const std::vector<std::set<IntPair>>& normalLandmarkGroups,
  const std::vector<std::set<IntPair>>& handleLandmarkGroups,
  VertexData<SurfacePoint>& AtoB,
  VertexData<SurfacePoint>& BtoA)
{
  SPDLOG_INFO(ANSI_BOLD ANSI_YELLOW "lift BEGIN" ANSI_RESET);

  lift_prepare(A_orig, B_orig, A_gmap, B_gmap, normalLandmarkGroups, handleLandmarkGroups);

  Pair<VertexData<SurfacePoint>&> vertexImage = { AtoB, BtoA };

  forEachModel2([&](Model& modelP, Model& modelQ) {
    SPDLOG_INFO("Computing vertex image from {} to {}", modelP.name, modelQ.name);

    // Stage I: Run path lifting on every vertex
    VertexData<std::vector<std::vector<SurfacePoint>>> pathQ(*modelP.orig->mesh);

    for (Vertex vP_orig : modelP.orig->mesh->vertices()) {
      Vertex vP_cut = util::convert_vertex(vP_orig, *modelP.gmap->mesh);

      std::vector<Halfedge> pathP;
      lift_main_root(modelP, modelQ, vP_cut, pathP, pathQ[vP_orig]);

      if (vP_orig.getIndex() && vP_orig.getIndex() % 1000 == 0)
        SPDLOG_INFO("  Processed {}", vP_orig);
    }

    // Stage II: Consolidate path lifting results based on conformal distortion
    SPDLOG_INFO("  Consolidating possible candidate solutions...");
    while (true) {
      // Pick a vertex to consolidate next which has maximum number of triangles each of which has two other vertices consolidated
      kt84::MaxSelector<Vertex> vNext;
      for (Vertex v : modelP.orig->mesh->vertices()) {
        // If this vertex has a single path lifting result, it is done
        if (pathQ[v].size() > 1) {
          // Count adjacent triangles whose two other vertices are already consolidated
          int n = 0;
          for (Halfedge he : v.outgoingHalfedges()) {
            Vertex v1 = he.tipVertex();
            Vertex v2 = he.next().tipVertex();
            if (pathQ[v1].size() == 1 && pathQ[v2].size() == 1)
              ++n;
          }
          vNext.update(n, v);
        }
      }

      // If all vertices are consolidated, we're done
      if (!vNext.count)
        break;

      // For the selected vertex, pick a solution that has minimal conformal distortion over all adjacent consolidated triangles
      kt84::MinSelector<int> bestIndex;
      for (size_t k = 0; k < pathQ[vNext.value].size(); ++k) {
        double totalDistortion = 0.0;

        for (Halfedge he : vNext.value.outgoingHalfedges()) {
          Vertex v1 = he.tipVertex();
          Vertex v2 = he.next().tipVertex();

          if (pathQ[v1].size() == 1 && pathQ[v2].size() == 1) {
            // Get triangle's corner positions in modelP.orig
            Vector3d p0 = kt84::vector_cast<3, Vector3d>(modelP.orig->geom->vertexPositions[vNext.value]);
            Vector3d p1 = kt84::vector_cast<3, Vector3d>(modelP.orig->geom->vertexPositions[v1]);
            Vector3d p2 = kt84::vector_cast<3, Vector3d>(modelP.orig->geom->vertexPositions[v2]);

            // Get triangle's corner positions when mapped onto modelQ.orig
            SurfacePoint sp0 = util::convert_sp(pathQ[vNext.value][k].back(), *modelQ.orig->mesh);
            SurfacePoint sp1 = util::convert_sp(pathQ[v1][0].back(), *modelQ.orig->mesh);
            SurfacePoint sp2 = util::convert_sp(pathQ[v2][0].back(), *modelQ.orig->mesh);

            Vector3d q0 = kt84::vector_cast<3, Vector3d>(sp0.interpolate(modelQ.orig->geom->vertexPositions));
            Vector3d q1 = kt84::vector_cast<3, Vector3d>(sp1.interpolate(modelQ.orig->geom->vertexPositions));
            Vector3d q2 = kt84::vector_cast<3, Vector3d>(sp2.interpolate(modelQ.orig->geom->vertexPositions));

            totalDistortion += util::compute_conformal_energy({p0, p1, p2}, {q0, q1, q2});
          }
        }

        bestIndex.update(totalDistortion, k);
      }

      pathQ[vNext.value] = { pathQ[vNext.value][bestIndex.value] };
    }

    // Return resulting vertex image
    modelP.select(vertexImage) = VertexData<SurfacePoint>(*modelP.orig->mesh);

    for (Vertex vP_orig : modelP.orig->mesh->vertices()) {
      SurfacePoint spQ_cut = pathQ[vP_orig][0].back();

      SurfacePoint spQ_orig = util::convert_sp(spQ_cut, *modelQ.orig->mesh);

      modelP.select(vertexImage)[vP_orig] = spQ_orig;
    }
  });

  SPDLOG_INFO(ANSI_BOLD ANSI_YELLOW "lift END" ANSI_RESET);
}

void lift_PtoQ(
  const MeshGeometry& A_orig, const MeshGeometry& B_orig,
  const MeshGeometry& A_gmap, const MeshGeometry& B_gmap,
  const std::vector<std::set<IntPair>>& normalLandmarkGroups,
  const std::vector<std::set<IntPair>>& handleLandmarkGroups,
  Vertex vP_orig,
  std::vector<Halfedge>& pathP_cut,
  std::vector<std::vector<SurfacePoint>>& pathQ_orig)
{
  lift_prepare(A_orig, B_orig, A_gmap, B_gmap, normalLandmarkGroups, handleLandmarkGroups);

  ASSERT_WITH_LOG(vP_orig.getMesh() == A_orig.mesh.get() || vP_orig.getMesh() == B_orig.mesh.get(), "vSource must come from either A_orig or B_orig");

  Model& modelP = vP_orig.getMesh() == A_orig.mesh.get() ? modelA : modelB;
  Model& modelQ = vP_orig.getMesh() == A_orig.mesh.get() ? modelB : modelA;

  Vertex vP_cut = util::convert_vertex(vP_orig, *modelP.gmap->mesh);

  Vertex vQ_cut = modelP.correspondingLandmark[vP_cut];

  if (vQ_cut == Vertex()) {
    std::vector<std::vector<SurfacePoint>> pathQ_cut;
    lift_main_root(modelP, modelQ, vP_cut, pathP_cut, pathQ_cut);

    // Convert pathQ_cut to pathQ_orig
    pathQ_orig.resize(pathQ_cut.size());
    for (size_t k = 0; k < pathQ_cut.size(); ++k) {
      pathQ_orig[k].resize(pathQ_cut[k].size());
      for (size_t i = 0; i < pathQ_cut[k].size(); ++i)
        pathQ_orig[k][i] = util::convert_sp(pathQ_cut[k][i], *modelQ.orig->mesh);
    }

  } else {
    pathP_cut = {};
    pathQ_orig = { { util::convert_sp(SurfacePoint(vQ_cut), *modelQ.orig->mesh) } };
  }

}

void lift_AtoB(
  const MeshGeometry& A_orig, const MeshGeometry& B_orig,
  const MeshGeometry& A_gmap, const MeshGeometry& B_gmap,
  const std::vector<std::set<IntPair>>& normalLandmarkGroups,
  const std::vector<std::set<IntPair>>& handleLandmarkGroups,
  Vertex vA_orig,
  std::vector<Halfedge>& pathA_cut,
  std::vector<std::vector<SurfacePoint>>& pathB_orig)
{
  ASSERT_WITH_LOG(vA_orig.getMesh() == A_orig.mesh.get(), "vA_orig must come from A_orig.mesh");

  lift_PtoQ(A_orig, B_orig, A_gmap, B_gmap, normalLandmarkGroups, handleLandmarkGroups, vA_orig, pathA_cut, pathB_orig);
}

void lift_BtoA(
  const MeshGeometry& A_orig, const MeshGeometry& B_orig,
  const MeshGeometry& A_gmap, const MeshGeometry& B_gmap,
  const std::vector<std::set<IntPair>>& normalLandmarkGroups,
  const std::vector<std::set<IntPair>>& handleLandmarkGroups,
  Vertex vB_orig,
  std::vector<Halfedge>& pathB_cut,
  std::vector<std::vector<SurfacePoint>>& pathA_orig)
{
  ASSERT_WITH_LOG(vB_orig.getMesh() == B_orig.mesh.get(), "vB_orig must come from B_orig.mesh");

  lift_PtoQ(A_orig, B_orig, A_gmap, B_gmap, normalLandmarkGroups, handleLandmarkGroups, vB_orig, pathB_cut, pathA_orig);
}
