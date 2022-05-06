#include "common.h"

void gmap(
  // input
  ManifoldSurfaceMesh& A_mesh,
  ManifoldSurfaceMesh& B_mesh,
  VertexPositionGeometry& A_geom_3D,
  VertexPositionGeometry& B_geom_3D,
  const std::vector<std::set<IntPair>>& normalLandmarkGroups,
  const std::vector<std::set<IntPair>>& handleLandmarkGroups,
  int num_iter,
  // output
  std::unique_ptr<VertexPositionGeometry>& A_geom_2D,
  std::unique_ptr<VertexPositionGeometry>& B_geom_2D,
  bool debugOutput = false);
