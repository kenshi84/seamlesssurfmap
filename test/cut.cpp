#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/mesh_graph_algorithms.h>

#include <set>
#include <map>

#include "util.h"

int main(int argc, char *argv[]) {
  /*
  +- Terminology -------------------------------------------------------------------------+
  |                                                                                       |
  | Handle landmarks:                                                                     |
  |   Specified on each handle. After cutting the high genus surface into genus zero,     |
  |   each cut loop corresponds to a pair of coincident boundary loops. Each handle       |
  |   landmark must be specified as a coincident vertex pair on such a pair of coincident |
  |   boundary loops.                                                                     |
  |                                                                                       |
  | Normal landmarks:                                                                     |
  |   Specified elsewhere individually.                                                   |
  |                                                                                       |
  +---------------------------------------------------------------------------------------+
  */

  // Load a manifold surface mesh from file
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;

#if 0
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/A-orig.obj");
  const std::set<int> landmarks = { 5105, 5433, 5805, 6722, 14143, 14839, 15941, 17597 };
#endif

#if 0
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/B-orig.obj");
  const std::set<int> landmarks = { 8713, 4828, 6294, 2669, 3288, 2977, 7543, 1643 };
#endif

#if 0
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/vase/A-genus0.obj");
  const std::set<int> landmarks = { 111, 1111, 1224, 1526, 2431, /* handle landmark pair */ 1268, 2518 };
#endif

#if 1
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/vase/B-genus0.obj");
  const std::set<int> landmarks = { 1899, 9, 2165, 409, 405, /* handle landmark pair */ 1014, 2538 };
#endif

#if 0
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pretzel-tripletorus/A-genus0.obj");
  const std::set<int> landmarks = { 201, 496, 952, 1463, 1507, 1, 847, 1332, 1986, 2006, 2026 };
#endif

#if 0
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pretzel-tripletorus/B-genus0.obj");
  const std::set<int> landmarks = { 3635, 1507, 2118, 4063, 1329, 3248, 2398, 1666, 4590, 4622, 4654 };
#endif

  std::unique_ptr<ManifoldSurfaceMesh> meshOrig;
  std::unique_ptr<VertexPositionGeometry> geometryOrig;
  util::copy_mesh(*mesh, *geometry, meshOrig, geometryOrig);

  // Assume no handle
  ASSERT(mesh->genus() == 0);
  ASSERT(mesh->nConnectedComponents() == 1);
  ASSERT(mesh->nBoundaryLoops() % 2 == 0);

  const int nHandles = mesh->nBoundaryLoops() / 2;

  // Establish boundary correspondence while checking consistency
  std::map<BoundaryLoop, BoundaryLoop> boundaryLoopMap;
  std::map<Vertex, Vertex> boundaryVertexMap;

  auto getDistance = [&](Vertex vA, Vertex vB) {
    Vector3 pA = geometry->vertexPositions[vA];
    Vector3 pB = geometry->vertexPositions[vB];
    return norm(pA - pB);
  };

  for (BoundaryLoop blA : mesh->boundaryLoops()) {
    for (BoundaryLoop blB : mesh->boundaryLoops()) {
      // Skip redundant tests
      if (blA >= blB) continue;

      // Early out
      if (blA.degree() != blB.degree()) continue;

      // Check distance to nearest vertex on blB
      kt84::MinSelector<Vertex> nearestVertex;
      for (Vertex vB : blB.adjacentVertices())
        nearestVertex.update(getDistance(blA.halfedge().vertex(), vB), vB);

      // Skip if the distance is not zero
      if (nearestVertex.score > 0) continue;

      // Matching loop found!
      Halfedge heA = blA.halfedge();
      Halfedge heB = nearestVertex.value.halfedge().twin().next();

      do {
        Vertex vA = heA.vertex();
        Vertex vB = heB.vertex();

        ASSERT(getDistance(vA, vB) == 0);

        boundaryVertexMap[vA] = vB;
        boundaryVertexMap[vB] = vA;

        heA = heA.next();
        heB = heB.prevOrbitVertex();

      } while (heA != blA.halfedge());

      boundaryLoopMap[blA] = blB;
      boundaryLoopMap[blB] = blA;

      break;
    }
  }

  ASSERT((int)boundaryLoopMap.size() == nHandles * 2);

  std::map<Vertex, size_t> normalLandmarks;   // Mapping from normal landmark vertex to normal landmark group ID
  std::set<Vertex> handleLandmarksSet;

  for (int i : landmarks) {
    ASSERT(i < mesh->nVertices());
    Vertex v = mesh->vertex(i);

    if (v.isBoundary())
      handleLandmarksSet.insert(v);
    else
      normalLandmarks.insert({v, normalLandmarks.size()});
  }
  ASSERT(normalLandmarks.size() >= 3);
  ASSERT((int)handleLandmarksSet.size() == nHandles * 2);

  // Assign group ID to handle landmarks
  std::map<Vertex, size_t> handleLandmarks;
  BoundaryLoopData<Vertex> handleLandmarkPerBL(*mesh);
  for (int i = 0; i < nHandles; ++i) {
    // Pop the first item in handleLandmarksSet
    Vertex v = *handleLandmarksSet.begin();
    handleLandmarksSet.erase(v);

    // The other vertex should also be found in the set
    Vertex vOther = boundaryVertexMap.at(v);
    ASSERT(handleLandmarksSet.count(vOther));
    handleLandmarksSet.erase(vOther);

    // Each boundary should have exactly one landmark specified
    BoundaryLoop bl = v.halfedge().twin().face().asBoundaryLoop();
    BoundaryLoop blOther = boundaryLoopMap.at(bl);
    ASSERT(handleLandmarkPerBL[bl] == Vertex());
    ASSERT(handleLandmarkPerBL[blOther] == Vertex());
    handleLandmarkPerBL[bl] = v;
    handleLandmarkPerBL[blOther] = vOther;

    handleLandmarks.insert({v, handleLandmarks.size()});
    handleLandmarks.insert({vOther, handleLandmarks.size()});
  }

  // Unique number of distinct landmarks (grouping duplicated vertices as one)
  size_t nNormalLandmarkGroups = normalLandmarks.size();
  size_t nHandleLandmarkGroups = handleLandmarks.size();

  // Pick the normal landmark pair whose shortest path is smallest
  kt84::MinSelector<std::vector<Halfedge>> shortestCutPath;
  for (auto v0 = normalLandmarks.begin(); v0 != normalLandmarks.end(); ++v0) {
    auto v1 = v0;
    for (++v1; v1 != normalLandmarks.end(); ++v1) {
      std::vector<Halfedge> path = shortestEdgePath(*geometry, v0->first, v1->first, true);
      ASSERT(!path.empty());

      double pathLength = util::get_path_length(*geometry, path);
      shortestCutPath.update(pathLength, path);
    }
  }
  SPDLOG_INFO("Initial cut path: from {} to {}",
    shortestCutPath.value.front().tailVertex(),
    shortestCutPath.value.back().tipVertex());

  // Cut along the path
  std::map<Vertex, Vertex> vertexMap;
  util::cut_along_path(*mesh, *geometry, shortestCutPath.value, vertexMap);
  util::restore_face_vertex_ordering(*geometryOrig, *geometry);

  SPDLOG_INFO("Saving initial result to cut-init.obj");
  writeSurfaceMesh(*mesh, *geometry, "cut-init.obj");

  // Classify normal landmarks into processed (on the boundary) and unprocessed (in the interior)
  auto classify_landmarks = [&](std::set<Vertex>& processed_landmarks, std::set<Vertex>& unprocessed_landmarks) {
    processed_landmarks = unprocessed_landmarks = {};
    for (auto [v, i] : normalLandmarks) {
      (v.isBoundary() ? processed_landmarks : unprocessed_landmarks).insert(v);
    }
  };

  // Repeat until unprocessed landmark set becomes empty
  int count = 0;
  while (true) {
    std::set<Vertex> processed_landmarks;
    std::set<Vertex> unprocessed_landmarks;
    classify_landmarks(processed_landmarks, unprocessed_landmarks);
    if (unprocessed_landmarks.empty()) break;

    // Find the {processed, unprocessed} landmark pair with smallest shortest path
    SPDLOG_INFO("Finding the normal landmark pair with smallest shortest path length");
    shortestCutPath = {};
    for (Vertex processed_landmark : processed_landmarks) {
      for (Vertex unprocessed_landmark : unprocessed_landmarks) {
        SPDLOG_INFO("Examining pair: {} & {}", processed_landmark, unprocessed_landmark);
        std::vector<Halfedge> path = shortestEdgePath(*geometry, processed_landmark, unprocessed_landmark, true);
        if (path.empty()) {
          SPDLOG_INFO("No path found between the two vertices, skipping");
          continue;
        }

        double pathLength = util::get_path_length(*geometry, path);
        SPDLOG_INFO("Path length: {}; current shortest: {}", pathLength, shortestCutPath.score);
        shortestCutPath.update(pathLength, path);
      }
    }
    ASSERT(shortestCutPath.count);
    Vertex vFront = shortestCutPath.value.front().tailVertex();
    Vertex vBack  = shortestCutPath.value.back ().tipVertex ();
    SPDLOG_INFO("Shortest cut path: from {} to {}", vFront, vBack);

    // Cut along the path
    util::cut_along_path(*mesh, *geometry, shortestCutPath.value, vertexMap);
    util::restore_face_vertex_ordering(*geometryOrig, *geometry);

    // Save interme
    SPDLOG_INFO("Saving intermediate result to cut-{}.obj", count);
    writeSurfaceMesh(*mesh, *geometry, "cut-" + std::to_string(count) + ".obj");
    ++count;

    // Register newly created vertex as landmark
    ASSERT(vertexMap.count(vFront));    // We know vFront was on the cut graph, so a new vertex will necessarily be created
    ASSERT(!vertexMap.count(vBack));    // We know vBack was in interior, so no new vertex is created
    normalLandmarks.insert({vertexMap.at(vFront), normalLandmarks.at(vFront)});
  }

  // Connect handle landmarks to the cut graph
  SPDLOG_INFO("Connecting handle landmarks to the cut graph by shortest path");
  handleLandmarksSet = {};
  for (auto [vHandle, i] : handleLandmarks) {
    handleLandmarksSet.insert(vHandle);
  }
  for (Vertex vHandle : handleLandmarksSet) {
    shortestCutPath = {};
    for (auto [vNormal, i] : normalLandmarks) {
      SPDLOG_INFO("Examining pair: {} & {}", vHandle, vNormal);
      std::vector<Halfedge> path = shortestEdgePath(*geometry, vHandle, vNormal, true);
      if (path.empty()) {
        SPDLOG_INFO("No path found between the two vertices, skipping");
        continue;
      }

      double pathLength = util::get_path_length(*geometry, path);
      SPDLOG_INFO("Path length: {}; current shortest: {}", pathLength, shortestCutPath.score);
      shortestCutPath.update(pathLength, path);
    }
    ASSERT(shortestCutPath.count);

    Vertex vFront = shortestCutPath.value.front().tailVertex();
    Vertex vBack  = shortestCutPath.value.back ().tipVertex ();
    SPDLOG_INFO("Shortest cut path: from {} to {}", vFront, vBack);

    // Cut along the path
    util::cut_along_path(*mesh, *geometry, shortestCutPath.value, vertexMap);
    util::restore_face_vertex_ordering(*geometryOrig, *geometry);

    // Save interme
    SPDLOG_INFO("Saving intermediate result to cut-{}.obj", count);
    writeSurfaceMesh(*mesh, *geometry, "cut-" + std::to_string(count) + ".obj");
    ++count;

    // Since we're merging two loops into one, new vertices are always created at both endpoints of the path; register them as landmarks
    ASSERT(vertexMap.count(vFront));
    ASSERT(vertexMap.count(vBack ));
    handleLandmarks.insert({vertexMap.at(vFront), handleLandmarks.at(vFront)});
    normalLandmarks.insert({vertexMap.at(vBack ), normalLandmarks.at(vBack )});
  }

  // Sort grouped normal/handle landmarks
  std::vector<std::set<Vertex>> normalLandmarkGroups(nNormalLandmarkGroups);
  std::vector<std::set<Vertex>> handleLandmarkGroups(nHandleLandmarkGroups);
  for (auto [v, i] : normalLandmarks) normalLandmarkGroups[i].insert(v);
  for (auto [v, i] : handleLandmarks) handleLandmarkGroups[i].insert(v);

  // Output landmark information as string
  std::ostringstream oss;

  // For normal landmarks
  std::string separatorOuter = "";
  for (size_t i = 0; i < nNormalLandmarkGroups; ++i) {
    oss << separatorOuter;

    std::string separatorInner = "";
    for (auto v : normalLandmarkGroups[i]) {
      oss << separatorInner << v.getIndex();
      separatorInner = ",";
    }

    separatorOuter = ";";
  }

  if (nHandleLandmarkGroups)
    oss << "#";

  // For handle landmarks
  separatorOuter = "";
  for (size_t i = 0; i < nHandleLandmarkGroups; ++i) {
    oss << separatorOuter;

    ASSERT(handleLandmarkGroups[i].size() == 2);
    Vertex v0 = *handleLandmarkGroups[i].begin();
    Vertex v1 = *++handleLandmarkGroups[i].begin();
    oss << v0.getIndex() << "," << v1.getIndex();

    separatorOuter = ";";
  }

  SPDLOG_INFO("Landmark information: {}", oss.str());

  SPDLOG_INFO("Done!");
}
