#include <geometrycentral/surface/barycentric_coordinate_helpers.h>

#include "refine.h"
#include "util.h"

namespace {

struct SplitEdgeInfo {
  size_t fIndex;
  int heIndexInF;

  SplitEdgeInfo() {}

  SplitEdgeInfo(Halfedge he) :
    fIndex(he.face().getIndex()),
    heIndexInF(halfedgeIndexInTriangle(he))
  {}
};

void splitEdge(MeshGeometry& mg, const SplitEdgeInfo& sei) {
  // Retrieve halfedge to be split
  Halfedge he = halfedgeInTriangleByIndex(mg.mesh->face(sei.fIndex), sei.heIndexInF);

  // Set he's orientation to canonical
  he = he.edge().halfedge();

  // Take edge endpoints for later use
  Vector3 p0 = mg.geom->vertexPositions[he.tailVertex()];
  Vector3 p1 = mg.geom->vertexPositions[he.tipVertex()];

  // If UVs are present, copy them as well
  Vector2 qL0, qL1, qL2;
  Vector2 qR0, qR1, qR2;
  if (mg.uv) {
    qL0 = (*mg.uv)[he.corner()];
    qL1 = (*mg.uv)[he.next().corner()];
    qL2 = (*mg.uv)[he.next().next().corner()];

    qR0 = (*mg.uv)[he.twin().corner()];
    qR1 = (*mg.uv)[he.twin().next().corner()];
    qR2 = (*mg.uv)[he.twin().next().next().corner()];
  }

  // Perform split, get pointer to newly inserted vertex
  Vertex vNew = mg.mesh->splitEdgeTriangular(he.edge()).vertex();

  // Set position of the new vertex to the midpoint of the previous edge endpoints
  mg.geom->vertexPositions[vNew] = 0.5 * (p0 + p1);

  if (mg.uv) {
    (*mg.uv)[he.corner()]               = (*mg.uv)[he.prevOrbitFace().twin().corner()]      = 0.5 * (qL0 + qL1);
    (*mg.uv)[he.twin().next().corner()] = (*mg.uv)[he.twin().next().twin().next().corner()] = 0.5 * (qR0 + qR1);
    (*mg.uv)[he.prevOrbitFace().corner()] = qL2;
    (*mg.uv)[he.prevOrbitFace().twin().prevOrbitFace().corner()] = qL0;
    (*mg.uv)[he.twin().next().twin().corner()] = qR2;
  }
}

}

void refine(
  const std::set<int>& landmarks_i,
  MeshGeometry& orig,
  MeshGeometry& genus0,
  MeshGeometry& disk)
{
  SPDLOG_INFO(ANSI_BOLD ANSI_YELLOW "refine BEGIN" ANSI_RESET);

  // Topology sanity check
  ASSERT(orig.mesh->nConnectedComponents() == 1);
  ASSERT(orig.mesh->nBoundaryLoops() == 0);

  const int nHandles = orig.mesh->genus();

  if (nHandles > 0) {
    ASSERT(genus0.mesh && genus0.geom);
    ASSERT(genus0.mesh->nConnectedComponents() == 1);
    ASSERT(genus0.mesh->nBoundaryLoops() == 2 * nHandles);
    ASSERT(genus0.mesh->genus() == 0);
    ASSERT(util::compare_orig_vs_cut(*orig.mesh, *orig.geom, *genus0.mesh, *genus0.geom));

  } else {
    ASSERT(!genus0.mesh && !genus0.geom);
  }

  ASSERT(disk.mesh->nConnectedComponents() == 1);
  ASSERT(disk.mesh->nBoundaryLoops() == 1);
  ASSERT(disk.mesh->genus() == 0);
  ASSERT(util::compare_orig_vs_cut(*orig.mesh, *orig.geom, *disk.mesh, *disk.geom));

  // Convert from int-based to element-based
  std::set<Vertex> landmarks;
  for (int landmark_i : landmarks_i) {
    ASSERT(landmark_i < disk.mesh->nVertices());
    landmarks.insert(disk.mesh->vertex(landmark_i));
  }

  // Ensure that no two landmarks are directly adjacent
  for (Vertex landmark : landmarks) {
    for (Vertex v : landmark.adjacentVertices())
      ASSERT_WITH_LOG(!landmarks.count(v), "Landmark {} and {} are directly adjacent", landmark, v);
  }

  // Split internal edges connecting boundaries
  size_t nSplits = 0;
  while (true) {
    bool found = false;
    for (Edge e : disk.mesh->edges()) {
      // Skip boundary edge
      if (e.isBoundary()) continue;

      // Skip if not both of its two vertices are on the boundary
      std::array<Vertex, 2> v = e.adjacentVertices();
      if (!v[0].isBoundary() || !v[1].isBoundary()) continue;

      // Skip if any of its opposite vertices is a landmark (this edge will be split in the later process)
      if (landmarks.count(e.halfedge().next().tipVertex()) || landmarks.count(e.halfedge().twin().next().tipVertex())) continue;

      // Found relevant edge to split!
      found = true;
      ++nSplits;

      SplitEdgeInfo sei(e.halfedge());

      for (MeshGeometry* mg : {&orig, &genus0, &disk}) {
        // Skip genus0 when orig genus is zero
        if (nHandles == 0 && mg == &genus0) continue;

        splitEdge(*mg, sei);
      }

      break;
    }
    if (!found) break;
  }

  SPDLOG_INFO("Split {} edges connecting boundaries", nSplits);

  // Collect all edges on each landmark's 1-ring rim for splitting
  std::vector<SplitEdgeInfo> splitEdgeInfos;
  std::set<Edge> processedEdge;
  for (Vertex landmark : landmarks) {
    for (Halfedge he : landmark.outgoingHalfedges()) {
      if (!he.isInterior()) continue;

      // Get halfedge opposite to landmark
      he = he.next();

      // This edge must be interior due to the earlier edge splitting process
      ASSERT(!he.edge().isBoundary());

      // If the edge was already processed, skip
      if (processedEdge.count(he.edge())) continue;
      processedEdge.insert(he.edge());

      splitEdgeInfos.emplace_back(he);
    }
  }

  // Perform refinement
  for (SplitEdgeInfo sei : splitEdgeInfos) {
    for (MeshGeometry* mg : {&orig, &genus0, &disk}) {
      // Skip genus0 when orig genus is zero
      if (nHandles == 0 && mg == &genus0) continue;

      splitEdge(*mg, sei);
    }
  }

  SPDLOG_INFO("Split {} edges around landmarks", splitEdgeInfos.size());

  // Make meshes compressed
  for (MeshGeometry* mg : {&orig, &genus0, &disk}) {
    // Skip genus0 when orig genus is zero
    if (nHandles == 0 && mg == &genus0) continue;

    mg->mesh->compress();
  }

  if (nHandles > 0) {
    ASSERT(util::compare_orig_vs_cut(*orig.mesh, *orig.geom, *genus0.mesh, *genus0.geom));
  }
  ASSERT(util::compare_orig_vs_cut(*orig.mesh, *orig.geom, *disk.mesh, *disk.geom));

  SPDLOG_INFO(ANSI_BOLD ANSI_YELLOW "refine END" ANSI_RESET);
}
