#include "gmap.h"

int main() {
  // Load a manifold surface mesh from file
  MeshGeometry A_disk;
  MeshGeometry B_disk;

#if 0
  std::tie(A_disk.mesh, A_disk.geom) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/A-disk.obj");
  std::tie(B_disk.mesh, B_disk.geom) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/B-disk.obj");
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
  std::tie(A_disk.mesh, A_disk.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/A-disk-refined.obj");
  std::tie(B_disk.mesh, B_disk.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/B-disk-refined.obj");
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

  std::unique_ptr<VertexPositionGeometry> A_geom_2D;
  std::unique_ptr<VertexPositionGeometry> B_geom_2D;

  gmap(
    *A_disk.mesh,
    *B_disk.mesh,
    *A_disk.geom,
    *B_disk.geom,
    normalLandmarkGroups,
    handleLandmarkGroups,
    10000,
    A_geom_2D,
    B_geom_2D,
    true);

  SPDLOG_INFO("Writing result to gmapA.obj / gmapB.obj");
  writeSurfaceMesh(*A_disk.mesh, *A_geom_2D, "gmapA.obj");
  writeSurfaceMesh(*B_disk.mesh, *B_geom_2D, "gmapB.obj");
}
