#include "cocut.h"

int main() {
  // Load a manifold surface mesh from file
  MeshGeometry A_genus0;
  MeshGeometry B_genus0;

#if 0
  std::tie(A_genus0.mesh, A_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/A-orig.obj");
  std::tie(B_genus0.mesh, B_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/bumpy-round-cube/B-orig.obj");
  const std::set<IntPair> landmarks = {
    {5105, 8713},
    {5433, 4828},
    {5805, 6294},
    {6722, 2669},
    {14143, 3288},
    {14839, 2977},
    {15941, 7543},
    {17597, 1643}
  };
#endif

#if 1
  std::tie(A_genus0.mesh, A_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/A-genus0.obj");
  std::tie(B_genus0.mesh, B_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/B-genus0.obj");
  const std::set<IntPair> landmarks = {
    { 111 , 1899 },
    { 1111, 9    },
    { 1224, 2165 },
    { 1526, 409  },
    { 2431, 405  },
    /* handle landmark pair */
    { 1268, 2537 },
    { 2518, 1014 }
  };
#endif

  MeshGeometry A_disk;
  MeshGeometry B_disk;
  std::vector<std::set<IntPair>> normalLandmarkGroups;
  std::vector<std::set<IntPair>> handleLandmarkGroups;

  cocut(A_genus0, B_genus0, landmarks, A_disk, B_disk, normalLandmarkGroups, handleLandmarkGroups);

  SPDLOG_INFO("Writing result to cocutA.obj / cocutB.obj");
  writeSurfaceMesh(*A_disk.mesh, *A_disk.geom, "cocutA.obj");
  writeSurfaceMesh(*B_disk.mesh, *B_disk.geom, "cocutB.obj");

  // Output landmark information as string
  std::ostringstream oss;

  // For normal landmarks
  std::string separatorOuter = "";
  for (size_t g = 0; g < normalLandmarkGroups.size(); ++g) {
    oss << separatorOuter;

    std::string separatorInner = "";
    for (auto i : normalLandmarkGroups[g]) {
      oss << separatorInner << std::to_string(i);
      separatorInner = ",";
    }

    separatorOuter = ";";
  }

  // For handle landmarks
  separatorOuter = "#";
  for (size_t g = 0; g < handleLandmarkGroups.size(); ++g) {
    oss << separatorOuter;

    ASSERT(handleLandmarkGroups[g].size() == 2);
    IntPair i0 = *handleLandmarkGroups[g].begin();
    IntPair i1 = *++handleLandmarkGroups[g].begin();
    oss << std::to_string(i0) << "," << std::to_string(i1);

    separatorOuter = ";";
  }

  SPDLOG_INFO("Landmark information: {}", oss.str());

  SPDLOG_INFO("Done!");
}
