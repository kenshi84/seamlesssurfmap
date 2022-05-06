#include "common.h"

// Splits the triangles adjacent to each landmark into two at the landmark
// Assumes no two landmarks are adjacent
void refine(
  const std::set<int>& landmarks,     // Vertex indices based on those in disk
  MeshGeometry& orig,
  MeshGeometry& genus0,
  MeshGeometry& disk);
