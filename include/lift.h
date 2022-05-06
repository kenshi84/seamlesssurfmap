#include "common.h"

void lift(
  const MeshGeometry& A_orig, const MeshGeometry& B_orig,
  const MeshGeometry& A_gmap, const MeshGeometry& B_gmap,
  const std::vector<std::set<IntPair>>& normalLandmarkGroups,
  const std::vector<std::set<IntPair>>& handleLandmarkGroups,
  VertexData<SurfacePoint>& AtoB,
  VertexData<SurfacePoint>& BtoA);

// Lifting of individual vertices:
// This may result in multiple possible solutions due to the ambiguity caused by double-winded intervals of the
// boundary loop. To resolve this ambiguity, use the other method that maps all vertices at once.

void lift_AtoB(
  const MeshGeometry& A_orig, const MeshGeometry& B_orig,
  const MeshGeometry& A_gmap, const MeshGeometry& B_gmap,
  const std::vector<std::set<IntPair>>& normalLandmarkGroups,
  const std::vector<std::set<IntPair>>& handleLandmarkGroups,

  Vertex vA_orig,                                         // From A_orig.mesh

  std::vector<Halfedge>& pathA_cut,                       // Path on A_gmap.mesh, from (suitable) landmark nearest to vA_cut to vA_cut,
                                                          // where vA_cut = convert_vertex(vA_orig, *A_gmap.mesh)

  std::vector<std::vector<SurfacePoint>>& pathB_orig);    // Path on B_orig.mesh, from the corresponding landmark to the mapped (face) point

void lift_BtoA(
  const MeshGeometry& A_orig, const MeshGeometry& B_orig,
  const MeshGeometry& A_gmap, const MeshGeometry& B_gmap,
  const std::vector<std::set<IntPair>>& normalLandmarkGroups,
  const std::vector<std::set<IntPair>>& handleLandmarkGroups,

  Vertex vB_orig,                                         // From B_orig.mesh

  std::vector<Halfedge>& pathB_cut,                       // Path on B_gmap.mesh, from (suitable) landmark nearest to vB_cut to vB_cut,
                                                          // where vB_cut = convert_vertex(vB_orig, *B_gmap.mesh)

  std::vector<std::vector<SurfacePoint>>& pathA_orig);    // Path on A_orig.mesh, from the corresponding landmark to the mapped (face) point
