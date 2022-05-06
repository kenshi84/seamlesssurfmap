#include "util.h"
#include "handlecut.h"

#include <kt84/ScopeExit.hh>

using namespace geometrycentral;
using namespace geometrycentral::surface;

bool handlecut(
  const MeshGeometry& input,
  int seedVertex,
  MeshGeometry& output,
  bool tunnelMode,
  int* correspondingVertex,
  handlecut_errorcode* ec)
{
  SPDLOG_INFO(ANSI_BOLD ANSI_YELLOW "handlecut BEGIN" ANSI_RESET);
  auto scopeExit = kt84::make_ScopeExit([](){
    SPDLOG_INFO(ANSI_BOLD ANSI_YELLOW "handlecut END" ANSI_RESET);
  });

  util::copy_mesh(input, output);

  ASSERT_WITH_LOG(seedVertex < output.mesh->nVertices(), "Seed vertex {} is no less than #vertices ({})", seedVertex, output.mesh->nVertices());
  const Vertex vSeed = output.mesh->vertex(seedVertex);
  SPDLOG_INFO("Seed vertex: {}", seedVertex);

  // Topological info before cutting
  const size_t nConnectedComponents_before = output.mesh->nConnectedComponents();
  const int genus_before = output.mesh->genus();
  SPDLOG_INFO("nConnectedComponents_before: {}", nConnectedComponents_before);
  SPDLOG_INFO("genus_before: {}", genus_before);
  ASSERT_WITH_LOG(genus_before != 0, "Mesh has no handle");

  // Some geometric info needed for determining handle loop
  output.geom->requireEdgeLengths();
  output.geom->requireVertexNormals();

  // Auxiliary data needed for Dijkstra
  VertexData<double> dist(*output.mesh, std::numeric_limits<double>::infinity());
  VertexData<Halfedge> prev(*output.mesh);

  // Trace shortest path from vSeed to v
  auto tracePath = [&](Vertex v) {
    std::vector<Halfedge> ret;
    while (v != vSeed) {
      ret.push_back(prev[v]);
      v = prev[v].tailVertex();
    }
    std::reverse(ret.begin(), ret.end());
    return ret;
  };

  std::set<Vertex> finalized = { vSeed };
  dist[vSeed] = 0;

  // Initialize candidates with direct neighbors
  std::set<Vertex> candidate;
  for (Halfedge heOut : vSeed.outgoingHalfedges()) {
    Vertex vNeighbor = heOut.tipVertex();
    candidate.insert(vNeighbor);
    dist[vNeighbor] = output.geom->edgeLengths[heOut.edge()];
    prev[vNeighbor] = heOut;
  }

  std::vector<Halfedge> cutPath;

  for (int iter = 0; !candidate.empty(); ++iter) {

    /*
    auto write_point_cloud = [&](std::set<Vertex>& vertexSet, const std::string& filename_base) {
      MatrixXd V(0, 3);
      for (Vertex v : vertexSet) {
        V.conservativeResize(V.rows() + 1, 3);
        V.row(V.rows() - 1) << kt84::vector_cast<3, Vector3d>(output.geom->vertexPositions[v]).transpose();
      }
      util::write_point_cloud(V, filename_base);
    };

    if (iter <= 200) {
      write_point_cloud(candidate, "candidate-" + std::to_string(iter));
      write_point_cloud(finalized, "finalized-" + std::to_string(iter));
    }
    */

    // Pick the current vertex in the candidate set with smallest distance
    kt84::MinSelector<Vertex> vCurr;
    for (Vertex v : candidate)
      vCurr.update(dist[v], v);

    // Move it from candidate to finalized
    candidate.erase(vCurr.value);
    finalized.insert(vCurr.value);

    for (Halfedge heOut : vCurr.value.outgoingHalfedges()) {
      Vertex vNeighbor = heOut.tipVertex();

      // Parent vertex should be skipped
      if (vNeighbor == prev[vCurr.value].tailVertex()) continue;

      // If vNeighbor is also finalized, check if the two shortest paths and heOut form a good loop
      if (finalized.count(vNeighbor)) {
        // Create tempoary cut path forming a loop
        std::vector<Halfedge> tempCutPath = tracePath(vCurr.value);
        tempCutPath.push_back(heOut);
        std::vector<Halfedge> otherPath = tracePath(vNeighbor);
        for (auto he = otherPath.rbegin(); he != otherPath.rend(); ++he) {
          tempCutPath.push_back(he->twin());
        }

        // Ensure the path passes each vertex exactly once
        std::set<Vertex> vertexSet;
        for (Halfedge he : tempCutPath) {
          vertexSet.insert(he.vertex());
        }
        if (vertexSet.size() == tempCutPath.size()) {
#if 0
          static int count = 0;
          std::vector<Vector3d> tempCutPathPoints(tempCutPath.size() + 1);
          for (size_t i = 0; i <= tempCutPath.size(); ++i)
            tempCutPathPoints[i] = kt84::vector_cast<3, Vector3d>(output.geom->vertexPositions[tempCutPath[i % tempCutPath.size()].vertex()]);
          util::write_polyline(tempCutPathPoints, "tempCutPath-" + std::to_string(count++));
#endif

          // Plausibility test based on estimated turning angle sum as well as sum of dot products with vertex normals
          double turningAngleSum = 0;
          Vector3 centroid = Vector3::zero();
          double lengthSum = 0;
          for (size_t i = 0; i < tempCutPath.size(); ++i) {
            // Two consecutive edge vectors and their lengths
            Halfedge he0 = tempCutPath[i];
            Halfedge he1 = tempCutPath[(i + 1) % tempCutPath.size()];

            Vector3 dir0 = output.geom->halfedgeVector(he0);
            Vector3 dir1 = output.geom->halfedgeVector(he1);

            double l0 = output.geom->edgeLengths[he0.edge()];
            double l1 = output.geom->edgeLengths[he1.edge()];

            // Position and normal at the center vertex
            Vertex v = he0.tipVertex();
            Vector3 n = output.geom->vertexNormals[v];
            Vector3 p = output.geom->vertexPositions[v];

            turningAngleSum += angleInPlane(dir0, dir1, n);
            centroid += 0.5 * (l0 + l1) * p;
            lengthSum += l0;
          }

          centroid /= lengthSum;

          // Deemed plausible when it's reasonably small (would be ~= 2pi for mostly planar loops)
          if (std::abs(turningAngleSum) < 0.5)  {
            // Compute winding number for handle-ness or tunnel-ness check
            int windingNumber = util::winding_number(*output.mesh, *output.geom, centroid);

            // We assume simple non self-intersecting watertight surface
            ASSERT(windingNumber == 0 || windingNumber == 1);

            if (!tunnelMode && windingNumber == 1 || tunnelMode && windingNumber == 0) {
              cutPath = tempCutPath;
              break;
            }
          }
        }
      }

      // Standard Dijkstra update process
      double distNew = dist[vCurr.value] + output.geom->edgeLengths[heOut.edge()];
      if (distNew < dist[vNeighbor]) {
        dist[vNeighbor] = distNew;
        prev[vNeighbor] = heOut;

        candidate.insert(vNeighbor);
      }
    }

    if (!cutPath.empty()) break;
  }

  if (cutPath.empty()) {
    SPDLOG_ERROR("Cut path not found");
    if (ec) *ec = HANDLECUT_PATH_NOTFOUND;
    return false;
  }

  // Cut along the path
  std::map<Vertex, Vertex> vertexMap;
  util::cut_along_path(*output.mesh, *output.geom, cutPath, vertexMap);
  util::restore_face_vertex_ordering(*input.geom, *output.geom);

  if (correspondingVertex) {
    *correspondingVertex = vertexMap.at(vSeed).getIndex();
  }

  // Topological info after cutting
  const size_t nConnectedComponents_after = output.mesh->nConnectedComponents();
  const int genus_after = output.mesh->genus();
  SPDLOG_INFO("nConnectedComponents_after: {}", nConnectedComponents_after);
  SPDLOG_INFO("genus_after: {}", genus_after);

  output.mesh->compress();

  // Check success
  if (nConnectedComponents_after != nConnectedComponents_before) {
    ASSERT(nConnectedComponents_after == nConnectedComponents_before + 1);      // The cut must have created a new connected component
    ASSERT(genus_after == genus_before);                                        // The genus should be the same after the cut

    SPDLOG_ERROR("The cut loop was invalid, creating a new connected component");
    if (ec) *ec = HANDLECUT_LOOP_INVALID;
    return false;
  }

  if (ec) *ec = HANDLECUT_SUCCESS;
  return true;
}
