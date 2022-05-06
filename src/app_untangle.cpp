#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <args/args.hxx>
#include <chrono>

#include "untangle.h"
#include "util.h"

int main(int argc, char *argv[]) {
  SPDLOG_INFO(ANSI_BOLD "Version: {}" ANSI_RESET, VERSIONTAG);

  // Setup args
  args::ArgumentParser parser("Mesh untangling app");
  args::HelpFlag help(parser, "help", "Display this help message", {'h', "help"});
  args::ValueFlag<double> argTheta(parser, "number", "Balance between angle preservation (0) and area preservation (1) (default: 1./128)", {"theta"}, 1./128);
  args::ValueFlag<int> argMaxIter(parser, "integer", "Maximum number of iterations (default: 10000)", {"maxiter"}, 10000);
  args::ValueFlag<double> argBFGSThreshold(parser, "number", "BFGS threshold (default: 1e-4)", {"bfgs-threshold"}, 1e-4);
  args::ValueFlag<int> argBFGSMaxIter(parser, "integer", "Maximum number of BFGS iterations (default: 30000)", {"maxiter"}, 30000);
  args::ValueFlag<double> argConvergenceThreshold(parser, "number", "Convergence threshold (default: 1e-5)", {"convergence-threshold"}, 1e-5);
  args::Flag argDebug(parser, "debug", "Print debug messages", {"debug"});
  args::Positional<std::string> argMesh3D(parser, "<mesh3d.obj>", "The input 3D mesh", args::Options::Required);
  args::Positional<std::string> argMesh2D(parser, "<mesh2d.obj>", "The input 2D mesh", args::Options::Required);
  args::Positional<std::string> argMeshOutput(parser, "<output.obj>", "The output mesh filename", args::Options::Required);

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::Error e) {
    SPDLOG_ERROR(e.what());
    std::cerr << parser;
    return 1;
  }

  MatrixXd V;
  MatrixXi F;
  if (!igl::read_triangle_mesh(args::get(argMesh3D), V, F)) {
    SPDLOG_ERROR("Failed to read {}", args::get(argMesh3D));
    return 1;
  }

  MatrixXd V_uv;
  MatrixXi F2;
  if (!igl::read_triangle_mesh(args::get(argMesh2D), V_uv, F2)) {
    SPDLOG_ERROR("Failed to read {}", args::get(argMesh2D));
    return 1;
  }

  ASSERT_WITH_LOG(V.rows() == V_uv.rows() && F.rows() == F2.rows(), "The 3D and 2D meshes must have the same connectivity");

  // Drop col(2) of V_uv
  VectorXd V_uv_col0 = V_uv.col(0);
  VectorXd V_uv_col1 = V_uv.col(1);
  V_uv.resize(V.rows(), 2);
  V_uv << V_uv_col0, V_uv_col1;

  UntangleParam param;
  param.theta = args::get(argTheta);
  param.maxiter = args::get(argMaxIter);
  param.bfgs_threshold = args::get(argBFGSThreshold);
  param.bfgs_maxiter = args::get(argBFGSMaxIter);
  param.debug = args::get(argDebug);
  param.converge_threshold = args::get(argConvergenceThreshold);

  if (!untangle(V, F, V_uv, &param)) {
    SPDLOG_ERROR("Failed to untangle!");
    return 1;
  }

  igl::write_triangle_mesh(args::get(argMeshOutput), util::append_zero_col(V_uv), F);
  return 0;
}
