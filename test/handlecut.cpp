#include "handlecut.h"

int main() {
  MeshGeometry input;

#if 0
  std::tie(input.mesh, input.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/A-orig.obj");
  const std::set<int> seedVertices = { 1268 };
#endif

#if 1
  std::tie(input.mesh, input.geom) = readManifoldSurfaceMesh(DATA_PATH "/vase/B-orig.obj");
  const std::set<int> seedVertices = { 1014 };
#endif


  MeshGeometry output;
  for (int i = 0; int seedVertex : seedVertices) {
    if (i++ > 0) {
      std::swap(input.mesh, output.mesh);
      std::swap(input.geom, output.geom);
    }
    ASSERT_WITH_LOG(handlecut(input, seedVertex, output), "Error in handlecut");
  }

  SPDLOG_INFO("Writing resulting mesh to handlecut.obj");
  writeSurfaceMesh(*output.mesh, *output.geom, "handlecut.obj");
}
