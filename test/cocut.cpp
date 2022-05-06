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

#if 0
  std::tie(A_genus0.mesh, A_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/cow-horse/A-orig.obj");
  std::tie(B_genus0.mesh, B_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/cow-horse/B-orig.obj");
  const std::set<IntPair> landmarks = {
    { 2, 2187 },
    { 44, 728 },
    { 65, 678 },
    { 257, 1759 },
    { 426, 449 },
    { 1021, 1867 },
    { 1537, 415 },
    { 1871, 1293 },
    { 1935, 1162 },
    { 2064, 1128 }
  };
#endif

#if 0
  std::tie(A_genus0.mesh, A_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/hands/A-orig.obj");
  std::tie(B_genus0.mesh, B_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/hands/B-orig.obj");
  const std::set<IntPair> landmarks = {
    { 85, 2922 },
    { 407, 23 },
    { 2703, 3494 },
    { 3192, 565 },
    { 3803, 4060 },
    { 3903, 2689 },
    { 4649, 1129 },
    { 6538, 351 },
    { 6722, 597 },
    { 7058, 1398 },
    { 8012, 4684 },
    { 8349, 3579 },
    { 8376, 433 },
    { 9037, 2057 },
    { 9491, 3419 },
    { 10759, 257 },
    { 10767, 3498 } 
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

#if 0
  std::tie(A_genus0.mesh, A_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pretzel-tripletorus/A-genus0.obj");
  std::tie(B_genus0.mesh, B_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pretzel-tripletorus/B-genus0.obj");
  const std::set<IntPair> landmarks = {
    { 356, 486 },
    { 358, 4200 },
    { 1246, 3452 },
    { 1380, 4322 },
    { 1416, 3612 },
    { 1420, 352 },
    { 1434, 4289 },
    { 1438, 3986 },
    { 1450, 3556 },
    { 1453, 885 },
    { 1489, 3433 },
    { 1526, 3689 },
    { 1593, 3456 },
    { 1685, 782 },
    { 1708, 1995 },
    { 1725, 3671 },
    { 1747, 1527 },
    { 1752, 3464 },
    { 547, 3247 },
    { 811, 2352 },
    { 1545, 1675 },
    { 1985, 4590 },
    { 2005, 4622 },
    { 2025, 4654 }
  };
#endif

#if 0
  std::tie(A_genus0.mesh, A_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pig-armadillo/A-orig.obj");
  std::tie(B_genus0.mesh, B_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/schmidt20/pig-armadillo/B-orig.obj");
  const std::set<IntPair> landmarks = {
    { 230, 2671 }, 
    { 291, 762 }, 
    { 389, 841 }, 
    { 814, 2077 }, 
    { 982, 1300 }, 
    { 1397, 2319 }, 
    { 2118, 633 }, 
    { 2474, 294 }, 
    { 2961, 1451 }, 
    { 3193, 343 }, 
    { 3534, 3377 }
  };
#endif

#if 0
  std::tie(A_genus0.mesh, A_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/../../cit-hot/data/schmidt20/duck-donut/A-genus0.obj");
  std::tie(B_genus0.mesh, B_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/../../cit-hot/data/schmidt20/duck-donut/B-genus0.obj");
  const std::set<IntPair> landmarks = {
    { 27, 5186 },
    { 216, 4035 },
    { 302, 320 },
    { 316, 3961 },
    { 459, 2552 },
    { 651, 423 },
    { 698, 516 },
    { 754, 1234 },
    { 918, 1187 },
    { 1013, 1214 },
    { 1067, 4970 },
    { 1095, 2748 },
    { 1170, 2368 },
    { 766, 3632 },
    { 1303, 5391 }
  };
#endif

#if 0
  std::tie(A_genus0.mesh, A_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/../../cit-hot/data/segeval/101-chair/105-102/var2/A-genus0.obj");
  std::tie(B_genus0.mesh, B_genus0.geom) = readManifoldSurfaceMesh(DATA_PATH "/../../cit-hot/data/segeval/101-chair/105-102/var2/B-genus0.obj");
  const std::set<IntPair> landmarks = {
    { 87, 409 },
    { 97, 2645 },
    { 112, 132 },
    { 155, 7027 },
    { 1007, 1366 },
    { 1125, 7324 },
    { 1272, 2012 },
    { 1282, 4978 },
    { 3108, 3925 },
    { 3194, 3918 },
    { 3203, 6372 },
    { 3215, 10369 },
    { 3217, 4581 },
    { 3674, 9278 },
    { 3688, 11475 },
    { 4513, 2523 },
    { 4750, 7179 },
    { 4925, 12925 },
    { 5113, 11365 },
    { 5354, 14739 },
    { 5409, 9513 },
    { 5492, 6214 },
    { 9233, 1647 },
    { 9234, 6543 },
    { 9246, 15498 },
    { 9260, 12858 },
    { 46, 216 },
    { 92, 1030 },
    { 102, 15798 },
    { 9295, 15760 },
    { 9330, 4043 },
    { 9365, 15833 }
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
