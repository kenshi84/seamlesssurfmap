#include "handlecut.h"

int main() {
  MeshGeometry input;
  std::tie(input.mesh, input.geom) = readManifoldSurfaceMesh(DATA_PATH "/pretzel.obj");

  const std::set<int> seedVertices_vaseA = { 1268 };
  const std::set<int> seedVertices_vaseB = { 1014 };
  const std::set<int> seedVertices_pretzel = { 784, 811, 1083 };
  const std::set<int> seedVertices_tripletorus = { 2392, 3234, 2891 };

  MeshGeometry output;
  for (int i = 0; int seedVertex : seedVertices_pretzel) {
    if (i++ > 0) {
      std::swap(input.mesh, output.mesh);
      std::swap(input.geom, output.geom);
    }
    ASSERT_WITH_LOG(handlecut(input, seedVertex, output), "Error in handlecut");
  }

  SPDLOG_INFO("Writing resulting mesh to handlecut.obj");
  writeSurfaceMesh(*output.mesh, *output.geom, "handlecut.obj");
}
