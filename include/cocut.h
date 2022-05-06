#include "common.h"

void cocut(
  const MeshGeometry& A_genus0,
  const MeshGeometry& B_genus0,
  const std::set<IntPair>& landmarks,
  MeshGeometry& A_disk,
  MeshGeometry& B_disk,
  std::vector<std::set<IntPair>>& normalLandmarkGroups,
  std::vector<std::set<IntPair>>& handleLandmarkGroups,
  bool debugOutput = false);
