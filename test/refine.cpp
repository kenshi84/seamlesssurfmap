#include "refine.h"

int main() {
  std::set<int> landmarks;
  MeshGeometry orig;
  MeshGeometry genus0;
  MeshGeometry disk;

#if 0
  std::tie(orig.mesh, orig.geom) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/A-orig.obj");
  std::tie(disk.mesh, disk.geom) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/A-disk.obj");
  landmarks = { 5105, 20247, 5433, 5805, 20013, 20105, 6722, 20199, 14143, 20059, 14839, 20151, 15941, 17597 };
#endif

#if 0
  std::tie(orig.mesh, orig.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/A-orig.obj");
  std::tie(genus0.mesh, genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/A-genus0.obj");
  std::tie(disk.mesh, disk.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/A-disk.obj");
  landmarks = { 111, 2565, 1111, 1224, 1526, 2532, 2547, 2431, 2618, 2641, 1268, 2594, 2518, 2619 };
#endif

#if 1
  std::tie(orig.mesh, orig.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/B-orig.obj");
  std::tie(genus0.mesh, genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/B-genus0.obj");
  std::tie(disk.mesh, disk.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/B-disk.obj");
  landmarks = { 1899, 2581, 9, 2165, 409, 2551, 2565, 405, 2629, 2656, 2537, 2607, 1014, 2630 };
#endif

#if 0
  std::tie(orig.mesh, orig.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pretzel-tripletorus/A-orig.obj");
  std::tie(genus0.mesh, genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pretzel-tripletorus/A-genus0.obj");
  std::tie(disk.mesh, disk.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pretzel-tripletorus/A-disk.obj");
  landmarks = { 356, 2160, 358, 2094, 2121, 1246, 2131, 1380, 2043, 1416, 2106, 1420, 2065, 2069, 1434, 1438, 1450, 2079, 1453, 2035, 2039, 2047, 1489, 2085, 2099, 1526, 2055, 1593, 2029, 1685, 2075, 2141, 1708, 2059, 2152, 1725, 1747, 2089, 1752, 2174, 547, 2111, 2005, 2153, 811, 2122, 2025, 2161, 1545, 2132, 1985, 2142 };
#endif

#if 0
  std::tie(orig.mesh, orig.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pig-armadillo/A-orig.obj");
  std::tie(disk.mesh, disk.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pig-armadillo/A-disk.obj");
  landmarks = { 230, 3758, 291, 3656, 389, 3610, 3633, 814, 3592, 982, 1397, 2118, 2474, 2961, 3679, 3701, 3801, 3193, 3534, 3723 };
#endif

  refine(landmarks, orig, genus0, disk);

  SPDLOG_INFO("Writing result to orig-refined.obj");
  writeSurfaceMesh(*orig.mesh, *orig.geom, "orig-refined.obj");

  if (genus0.mesh) {
    SPDLOG_INFO("Writing result to genus0-refined.obj");
    writeSurfaceMesh(*genus0.mesh, *genus0.geom, "genus0-refined.obj");
  }

  SPDLOG_INFO("Writing result to disk-refined.obj");
  writeSurfaceMesh(*disk.mesh, *disk.geom, "disk-refined.obj");
}
