#include <geometrycentral/surface/mesh_graph_algorithms.h>

#include "cocut.h"
#include "util.h"

namespace {

using Path = std::vector<Halfedge>;
using VertexMap = std::map<Vertex, Vertex>;
using BoundaryLoopMap = std::map<BoundaryLoop, BoundaryLoop>;

struct PathPair : public Pair<Path> {
  VertexPair vFront() const { return {A.front().tailVertex(), B.front().tailVertex()}; }
  VertexPair vBack() const { return {A.back().tipVertex(), B.back().tipVertex()}; }
};

struct VertexMapPair : public Pair<VertexMap> {
  VertexPair at(const VertexPair& v) const { return { A.at(v.A), B.at(v.B) }; }
  size_t count(const VertexPair& v) const {
    size_t cA = A.count(v.A);
    size_t cB = B.count(v.B);
    ASSERT(cA == cB);
    return cA;
  }
};

using BoundaryLoopMapPair = Pair<BoundaryLoopMap>;

template <typename T>
using BoundaryLoopDataPair = Pair<BoundaryLoopData<T>>;

struct Model : public ModelBase {
  const MeshGeometry* genus0;
  MeshGeometry* disk;
};

}

void cocut(
  const MeshGeometry& A_genus0,
  const MeshGeometry& B_genus0,
  const std::set<IntPair>& landmarks,
  MeshGeometry& A_disk,
  MeshGeometry& B_disk,
  std::vector<std::set<IntPair>>& normalLandmarkGroups,
  std::vector<std::set<IntPair>>& handleLandmarkGroups,
  bool debugOutput)
{
  SPDLOG_INFO(ANSI_BOLD ANSI_YELLOW "cocut BEGIN" ANSI_RESET);

  Model modelA = { "A", &A_genus0, &A_disk };
  Model modelB = { "B", &B_genus0, &B_disk };

  auto forEachModel = [&](std::function<void(Model& model)> f) {
    f(modelA);
    f(modelB);
  };

  forEachModel([&](Model& model) {
    // Copy input genus0 mesh into workig mesh which will eventually become disk topology
    util::copy_mesh(*model.genus0, *model.disk);

    ASSERT_WITH_LOG(model.disk->mesh->genus() == 0, "{}'s genus must be zero, but not ({})", model.name, model.disk->mesh->genus());
    ASSERT_WITH_LOG(model.disk->mesh->nConnectedComponents() == 1, "{}'s #connectedComponents must be one, but not ({})", model.name, model.disk->mesh->nConnectedComponents());
    ASSERT_WITH_LOG(model.disk->mesh->nBoundaryLoops() % 2 == 0, "{}'s #boundaryLoops must be even, but not ({})", model.name, model.disk->mesh->nBoundaryLoops());
  });
  ASSERT_WITH_LOG(modelA.disk->mesh->nBoundaryLoops() == modelB.disk->mesh->nBoundaryLoops(),
    "A and B have different number of boundaryLoops ({} vs {})",
    modelA.disk->mesh->nBoundaryLoops(),
    modelB.disk->mesh->nBoundaryLoops());

  const int nHandles = modelA.disk->mesh->nBoundaryLoops() / 2;

  BoundaryLoopMapPair boundaryLoopMap;
  VertexMapPair boundaryVertexMap;

  // Establish boundary correspondence while checking consistency
  forEachModel([&](Model& model) {
    auto getDistance = [&](Vertex v0, Vertex v1) {
      Vector3 p0 = model.disk->geom->vertexPositions[v0];
      Vector3 p1 = model.disk->geom->vertexPositions[v1];
      return norm(p0 - p1);
    };

    for (BoundaryLoop bl0 : model.disk->mesh->boundaryLoops()) {
      for (BoundaryLoop bl1 : model.disk->mesh->boundaryLoops()) {
        // Skip redundant tests
        if (bl0 >= bl1) continue;

        // Early out
        if (bl0.degree() != bl1.degree()) continue;

        // Check distance to nearest vertex on bl1
        kt84::MinSelector<Vertex> nearestVertex;
        for (Vertex v1 : bl1.adjacentVertices())
          nearestVertex.update(getDistance(bl0.halfedge().vertex(), v1), v1);

        // Skip if the distance is not zero
        if (nearestVertex.score > 0) continue;

        // Matching loop found!
        Halfedge he0 = bl0.halfedge();
        Halfedge he1 = nearestVertex.value.halfedge().twin().next();

        do {
          Vertex v0 = he0.vertex();
          Vertex v1 = he1.vertex();

          ASSERT(getDistance(v0, v1) == 0);

          model.select(boundaryVertexMap)[v0] = v1;
          model.select(boundaryVertexMap)[v1] = v0;

          he0 = he0.next();
          he1 = he1.prevOrbitVertex();

        } while (he0 != bl0.halfedge());

        model.select(boundaryLoopMap)[bl0] = bl1;
        model.select(boundaryLoopMap)[bl1] = bl0;

        break;
      }
    }

    ASSERT((int)model.select(boundaryLoopMap).size() == nHandles * 2);
  });

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

  std::map<VertexPair, size_t> normalLandmarks;   // Mapping from normal landmark vertex to normal landmark group ID
  std::set<VertexPair> handleLandmarksSet;

  for (auto [iA, iB] : landmarks) {
    ASSERT_WITH_LOG(iA < modelA.disk->mesh->nVertices(), "Specified vertex index {} is no less than A's nVertices {}", iA, modelA.disk->mesh->nVertices());
    ASSERT_WITH_LOG(iB < modelB.disk->mesh->nVertices(), "Specified vertex index {} is no less than B's nVertices {}", iB, modelB.disk->mesh->nVertices());
    VertexPair v = {
      modelA.disk->mesh->vertex(iA),
      modelB.disk->mesh->vertex(iB)
    };

    if (v.isBoundary())
      handleLandmarksSet.insert(v);
    else
      normalLandmarks.insert({v, normalLandmarks.size()});
  }
  ASSERT_WITH_LOG(normalLandmarks.size() >= 2, "At least 2 normal landmarks should be specified");
  ASSERT_WITH_LOG((int)handleLandmarksSet.size() == nHandles * 2, "Wrong number of handle landmarks {} (should be {})", handleLandmarksSet.size(), 2 * nHandles);

  // Assign group ID to handle landmarks
  std::map<VertexPair, size_t> handleLandmarks;
  BoundaryLoopDataPair<Vertex> handleLandmarkPerBL = {
    { *modelA.disk->mesh },
    { *modelB.disk->mesh }
  };
  for (int i = 0; i < nHandles; ++i) {
    // Pop the first item in handleLandmarksSet
    VertexPair v = *handleLandmarksSet.begin();
    handleLandmarksSet.erase(v);

    // The other vertex should also be found in the set
    VertexPair vOther = boundaryVertexMap.at(v);
    ASSERT_WITH_LOG(handleLandmarksSet.count(vOther), "Handle landmarks should be specified in pairs, but {}'s counterpart {} is not specified as landmark", v, vOther);
    handleLandmarksSet.erase(vOther);

    // Each boundary should have exactly one landmark specified
    forEachModel([&](Model& model) {
      BoundaryLoop bl = model.select(v).halfedge().twin().face().asBoundaryLoop();
      BoundaryLoop blOther = model.select(boundaryLoopMap).at(bl);
      ASSERT_WITH_LOG(model.select(handleLandmarkPerBL)[bl] == Vertex(), "{}'s boundaryLoop {} is specified multiple landmarks", model.name, bl.getIndex());
      ASSERT_WITH_LOG(model.select(handleLandmarkPerBL)[blOther] == Vertex(), "{}'s boundaryLoop {} is specified multiple landmarks", model.name, blOther.getIndex());
      model.select(handleLandmarkPerBL)[bl] = model.select(v);
      model.select(handleLandmarkPerBL)[blOther] = model.select(vOther);
    });

    handleLandmarks.insert({v, handleLandmarks.size()});
    handleLandmarks.insert({vOther, handleLandmarks.size()});
  }

  // Unique number of distinct landmarks (grouping duplicated vertices as one)
  size_t nNormalLandmarkGroups = normalLandmarks.size();
  size_t nHandleLandmarkGroups = handleLandmarks.size();

  // Pick the normal landmark pair whose shortest path is smallest
  kt84::MinSelector<PathPair> shortestCutPath;
  for (auto v0 = normalLandmarks.begin(); v0 != normalLandmarks.end(); ++v0) {
    auto v1 = v0;
    for (++v1; v1 != normalLandmarks.end(); ++v1) {
      PathPair path;
      double pathLength = 0;

      forEachModel([&](Model& model) {
        model.select(path) = shortestEdgePath(*model.disk->geom, model.select(v0->first), model.select(v1->first), true);
        ASSERT_WITH_LOG(!model.select(path).empty(), "Initial cut path between normal landmarks must be always found");

        pathLength += util::get_path_length(*model.disk->geom, model.select(path));
      });

      shortestCutPath.update(pathLength, path);
    }
  }

  // Cut along the path
  forEachModel([&](Model& model) {
    SPDLOG_INFO("Initial cut path: from {} to {} on {}",
      model.select(shortestCutPath.value.vFront()),
      model.select(shortestCutPath.value.vBack()),
      model.name);

    VertexMap vertexMap;
    util::cut_along_path(*model.disk->mesh, *model.disk->geom, model.select(shortestCutPath.value), vertexMap);
    util::restore_face_vertex_ordering(*model.genus0->geom, *model.disk->geom);

    if (debugOutput) {
      SPDLOG_INFO("Writing to cocut{}-0.obj", model.name);
      writeSurfaceMesh(*model.disk->mesh, *model.disk->geom, "cocut" + model.name + "-0.obj");
    }
  });

  size_t count = 1;

  // Repeat until unprocessed landmark set becomes empty
  while (true) {
    // Classify normal landmarks into processed (on the boundary) and unprocessed (in the interior)
    std::set<VertexPair> processed_landmarks;
    std::set<VertexPair> unprocessed_landmarks;
    for (auto [v, i] : normalLandmarks) {
      (v.isBoundary() ? processed_landmarks : unprocessed_landmarks).insert(v);
    }
    if (unprocessed_landmarks.empty()) break;

    // Find the {processed, unprocessed} landmark pair with smallest shortest path
    SPDLOG_INFO("Finding the normal landmark pair with smallest shortest path length");
    shortestCutPath = {};
    for (VertexPair processed_landmark : processed_landmarks) {
      for (VertexPair unprocessed_landmark : unprocessed_landmarks) {
        PathPair path;
        double pathLength = 0;

        forEachModel([&](Model& model) {
          if (pathLength == std::numeric_limits<double>::infinity()) return; // We already know that there's no shortest path on modelA, skip

          SPDLOG_INFO("Examining pair: {} & {} on {}",
            model.select(processed_landmark),
            model.select(unprocessed_landmark),
            model.name);

          model.select(path) = shortestEdgePath(*model.disk->geom,
            model.select(processed_landmark),
            model.select(unprocessed_landmark),
            true);

          if (model.select(path).empty()) {
            SPDLOG_INFO("No path found between the two vertices, skipping");
            pathLength = std::numeric_limits<double>::infinity();
            return;
          }

          pathLength += util::get_path_length(*model.disk->geom, model.select(path));
        });

        SPDLOG_INFO("Path length: {}; current shortest: {}", pathLength, shortestCutPath.score);
        shortestCutPath.update(pathLength, path);
      }
    }
    ASSERT(shortestCutPath.count);

    VertexPair vFront = shortestCutPath.value.vFront();
    VertexPair vBack = shortestCutPath.value.vBack();

    VertexMapPair vertexMap;

    forEachModel([&](Model& model) {
      SPDLOG_INFO("Shortest cut path: from {} to {} on {}",
        model.select(vFront),
        model.select(vBack),
        model.name);

      // Cut along the path
      util::cut_along_path(*model.disk->mesh, *model.disk->geom, model.select(shortestCutPath.value), model.select(vertexMap));
      util::restore_face_vertex_ordering(*model.genus0->geom, *model.disk->geom);

      if (debugOutput) {
        SPDLOG_INFO("Writing to cocut{}-{}.obj", model.name, count);
        writeSurfaceMesh(*model.disk->mesh, *model.disk->geom, "cocut" + model.name + "-" + std::to_string(count) + ".obj");
      }
    });

    ++count;

    // Register newly created vertex as landmark
    ASSERT(vertexMap.count(vFront));    // We know vFront was on the cut graph, so a new vertex will necessarily be created
    ASSERT(!vertexMap.count(vBack));    // We know vBack was in interior, so no new vertex is created
    normalLandmarks.insert({vertexMap.at(vFront), normalLandmarks.at(vFront)});
  }

  // Connect handle landmarks to the cut graph
  if (nHandleLandmarkGroups) {
    SPDLOG_INFO("Connecting handle landmarks to the cut graph by shortest path");
  }
  ASSERT(handleLandmarksSet.empty());
  for (auto [vHandle, i] : handleLandmarks) {
    handleLandmarksSet.insert(vHandle);
  }
  for (VertexPair vHandle : handleLandmarksSet) {
    shortestCutPath = {};
    for (auto [vNormal, i] : normalLandmarks) {
      PathPair path;
      double pathLength = 0;

      forEachModel([&, vNormal = vNormal /* https://stackoverflow.com/questions/50799719 */ ](Model& model) {
        if (pathLength == std::numeric_limits<double>::infinity()) return; // We already know that there's no shortest path on modelA, skip

        SPDLOG_INFO("Examining pair: {} & {} on {}",
          model.select(vHandle),
          model.select(vNormal),
          model.name);

        model.select(path) = shortestEdgePath(*model.disk->geom,
          model.select(vHandle),
          model.select(vNormal),
          true);

        if (model.select(path).empty()) {
          SPDLOG_INFO("No path found between the two vertices, skipping");
          pathLength = std::numeric_limits<double>::infinity();
          return;
        }

        pathLength += util::get_path_length(*model.disk->geom, model.select(path));
      });

      SPDLOG_INFO("Path length: {}; current shortest: {}", pathLength, shortestCutPath.score);
      shortestCutPath.update(pathLength, path);
    }
    ASSERT(shortestCutPath.count);

    VertexPair vFront = shortestCutPath.value.vFront();
    VertexPair vBack = shortestCutPath.value.vBack();

    VertexMapPair vertexMap;

    forEachModel([&](Model& model) {
      SPDLOG_INFO("Shortest cut path: from {} to {} on {}",
        model.select(vFront),
        model.select(vBack),
        model.name);

      // Cut along the path
      util::cut_along_path(*model.disk->mesh, *model.disk->geom, model.select(shortestCutPath.value), model.select(vertexMap));
      util::restore_face_vertex_ordering(*model.genus0->geom, *model.disk->geom);

      if (debugOutput) {
        SPDLOG_INFO("Writing to cocut{}-{}.obj", model.name, count);
        writeSurfaceMesh(*model.disk->mesh, *model.disk->geom, "cocut" + model.name + "-" + std::to_string(count) + ".obj");
      }
    });

    ++count;

    // Since we're merging two loops into one, new vertices are always created at both endpoints of the path; register them as landmarks
    ASSERT(vertexMap.count(vFront));
    ASSERT(vertexMap.count(vBack));
    handleLandmarks.insert({vertexMap.at(vFront), handleLandmarks.at(vFront)});
    normalLandmarks.insert({vertexMap.at(vBack), normalLandmarks.at(vBack)});
  }

  // Sort grouped normal/handle landmarks
  normalLandmarkGroups = std::vector<std::set<IntPair>>(nNormalLandmarkGroups);
  handleLandmarkGroups = std::vector<std::set<IntPair>>(nHandleLandmarkGroups);
  for (auto [v, i] : normalLandmarks) normalLandmarkGroups[i].insert(v.getIndex());
  for (auto [v, i] : handleLandmarks) handleLandmarkGroups[i].insert(v.getIndex());

  forEachModel([&](Model& model) {
    model.disk->mesh->compress();
  });

  SPDLOG_INFO(ANSI_BOLD ANSI_YELLOW "cocut END" ANSI_RESET);
}
