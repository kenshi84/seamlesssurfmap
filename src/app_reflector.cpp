#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>

#include "common.h"

int main(int argc, char *argv[]) {
  SPDLOG_INFO(ANSI_BOLD "Version: {}" ANSI_RESET, VERSIONTAG);

  MatrixXd V;
  MatrixXi F;
  igl::read_triangle_mesh(argv[1], V, F);

  V.col(0) *= -1;

  VectorXi col0 = F.col(0);
  VectorXi col2 = F.col(2);
  F.col(0) = col2;
  F.col(2) = col0;

  igl::write_triangle_mesh(argv[1] + std::string("-reflected.obj"), V, F);
}
