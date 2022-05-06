#include "lift.h"
#include "util.h"

int main() {
  MeshGeometry A_orig, B_orig;
  MeshGeometry A_gmap, B_gmap;

#if 0
  std::tie(A_orig.mesh, A_orig.geom) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/A-orig.obj");
  std::tie(B_orig.mesh, B_orig.geom) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/B-orig.obj");
  std::tie(A_gmap.mesh, A_gmap.geom) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/A-gmap.obj");
  std::tie(B_gmap.mesh, B_gmap.geom) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/B-gmap.obj");
  const std::vector<std::set<IntPair>> normalLandmarkGroups = {
    { { 5105, 8713 }, { 20247, 9202 } },
    { { 5433, 4828 } },
    { { 5805, 6294 }, { 20013, 9048 }, { 20105, 9108 } },
    { { 6722, 2669 }, { 20199, 9170 } },
    { { 14143, 3288 }, { 20059, 9075 } },
    { { 14839, 2977 }, { 20151, 9138 } },
    { { 15941, 7543 } },
    { { 17597, 1643 } }
  };
  const std::vector<std::set<IntPair>> handleLandmarkGroups = {};
#endif

#if 1
  std::tie(A_orig.mesh, A_orig.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/A-orig-refined.obj");
  std::tie(B_orig.mesh, B_orig.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/B-orig-refined.obj");
  std::tie(A_gmap.mesh, A_gmap.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/A-gmap.obj");
  std::tie(B_gmap.mesh, B_gmap.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/B-gmap.obj");
  const std::vector<std::set<IntPair>> normalLandmarkGroups = {
    { { 111, 1899 }, { 2565, 2581 } },
    { { 1111, 9 } },
    { { 1224, 2165 } },
    { { 1526, 409 }, { 2532, 2551 }, { 2547, 2565 } },
    { { 2431, 405 }, { 2618, 2629 }, { 2641, 2656 } }
  };
  const std::vector<std::set<IntPair>> handleLandmarkGroups = {
    { { 1268, 2537}, { 2594, 2607 } },
    { { 2518, 1014}, { 2619, 2630 } }
  };
#endif

#if 1
  // Test mapping of all vertices

  VertexData<SurfacePoint> AtoB;
  VertexData<SurfacePoint> BtoA;

  lift(A_orig, B_orig,
       A_gmap, B_gmap,
       normalLandmarkGroups,
       handleLandmarkGroups,
       AtoB, BtoA);

  SPDLOG_INFO("Writing result to AtoB.txt / BtoA.txt");

  std::ofstream fout("AtoB.txt");
  for (Vertex vA : A_orig.mesh->vertices()) {
    SurfacePoint spB = AtoB[vA];
    for (int i = 0; Vertex vB : spB.face.adjacentVertices()) {
      if (spB.faceCoords[i] > 0.0) {
        fout << vA.getIndex() << " " << vB.getIndex() << " " << spB.faceCoords[i] << "\n";
      }
      ++i;
    }
  }

  fout = std::ofstream("BtoA.txt");
  for (Vertex vB : B_orig.mesh->vertices()) {
    SurfacePoint spA = BtoA[vB];
    for (int i = 0; Vertex vA : spA.face.adjacentVertices()) {
      if (spA.faceCoords[i] > 0.0) {
        fout << vB.getIndex() << " " << vA.getIndex() << " " << spA.faceCoords[i] << "\n";
      }
      ++i;
    }
  }

  SPDLOG_INFO("Writing to AonB.obj / BonA.obj");

  VertexPositionGeometry AonB_geom(*A_orig.mesh);
  VertexPositionGeometry BonA_geom(*B_orig.mesh);
  for (Vertex v : A_orig.mesh->vertices()) AonB_geom.inputVertexPositions[v] = AtoB[v].interpolate(B_orig.geom->vertexPositions);
  for (Vertex v : B_orig.mesh->vertices()) BonA_geom.inputVertexPositions[v] = BtoA[v].interpolate(A_orig.geom->vertexPositions);
  writeSurfaceMesh(*A_orig.mesh, AonB_geom, "AonB.obj");
  writeSurfaceMesh(*B_orig.mesh, BonA_geom, "BonA.obj");

#endif

#if 0
  // Test mapping of single vertex (A -> B)

  Vertex vA_orig = A_orig.mesh->vertex(2203);

  std::vector<Halfedge> pathA_cut;
  std::vector<std::vector<SurfacePoint>> pathB_orig;

  lift_AtoB(
    A_orig, B_orig,
    A_gmap, B_gmap,
    normalLandmarkGroups,
    handleLandmarkGroups,
    vA_orig, pathA_cut, pathB_orig);

  // Write resulting paths as polylines
  MatrixXd V_orig;
  MatrixXd V_cut;

  std::ofstream pathA_cut_vertices("path_A_cut_vertices.txt");

  V_orig.resize(0, 3);
  V_cut.resize(0, 3);
  for (size_t i = 0; i < pathA_cut.size(); ++i) {
    Halfedge heA_cut = pathA_cut[i];

    for (int j = 0; (i == 0 && j < 2) || j < 1; ++j) {
      Vertex vA_cut = i == 0 && j == 0 ? heA_cut.tailVertex() : heA_cut.tipVertex();
      Vertex vA_orig = util::convert_vertex(vA_cut, *A_orig.mesh);

      V_orig.conservativeResize(V_orig.rows() + 1, 3);
      V_cut.conservativeResize(V_cut.rows() + 1, 3);
      V_orig.bottomRows(1) = kt84::vector_cast<3, RowVector3d>(A_orig.geom->vertexPositions[vA_orig]);
      V_cut.bottomRows(1) = kt84::vector_cast<3, RowVector3d>(A_gmap.geom->vertexPositions[vA_cut]);

      pathA_cut_vertices << vA_cut.getIndex() << std::endl;
    }
  }
  util::write_polyline(V_orig, "path_A_orig");
  util::write_polyline(V_cut, "path_A_gmap");

  for (size_t k = 0; k < pathB_orig.size(); ++k) {
    V_orig.resize(0, 3);
    for (SurfacePoint spB_orig : pathB_orig[k]) {
      V_orig.conservativeResize(V_orig.rows() + 1, 3);
      V_orig.bottomRows(1) = kt84::vector_cast<3, RowVector3d>(spB_orig.interpolate(B_orig.geom->vertexPositions));
    }
    util::write_polyline(V_orig, "path_B_orig-" + std::to_string(k));
  }
#endif

#if 0
  // Test mapping of single vertex (B -> A)

  Vertex vB_orig = B_orig.mesh->vertex(4574);

  std::vector<Halfedge> pathB_cut;
  std::vector<SurfacePoint> pathA_orig;

  lift_BtoA(
    A_orig, B_orig,
    A_gmap, B_gmap,
    normalLandmarkGroups,
    handleLandmarkGroups,
    vB_orig, pathB_cut, pathA_orig);

  // Write resulting paths as polylines
  MatrixXd V_orig;
  MatrixXd V_cut;

  std::ofstream pathB_cut_vertices("path_B_cut_vertices.txt");

  V_orig.resize(0, 3);
  V_cut.resize(0, 3);
  for (size_t i = 0; i < pathB_cut.size(); ++i) {
    Halfedge heB_cut = pathB_cut[i];

    for (int j = 0; (i == 0 && j < 2) || j < 1; ++j) {
      Vertex vB_cut = i == 0 && j == 0 ? heB_cut.tailVertex() : heB_cut.tipVertex();
      Vertex vB_orig = util::convert_vertex(vB_cut, *B_orig.mesh);

      V_orig.conservativeResize(V_orig.rows() + 1, 3);
      V_cut.conservativeResize(V_cut.rows() + 1, 3);
      V_orig.bottomRows(1) = kt84::vector_cast<3, RowVector3d>(B_orig.geom->vertexPositions[vB_orig]);
      V_cut.bottomRows(1) = kt84::vector_cast<3, RowVector3d>(B_gmap.geom->vertexPositions[vB_cut]);

      pathB_cut_vertices << vB_cut.getIndex() << std::endl;
    }
  }
  util::write_polyline(V_orig, "path_B_orig");
  util::write_polyline(V_cut, "path_B_gmap");

  V_orig.resize(0, 3);
  for (SurfacePoint spA_orig : pathA_orig) {
    V_orig.conservativeResize(V_orig.rows() + 1, 3);
    V_orig.bottomRows(1) = kt84::vector_cast<3, RowVector3d>(spA_orig.interpolate(A_orig.geom->vertexPositions));
  }
  util::write_polyline(V_orig, "path_A_orig");
#endif

}
