#include <args/args.hxx>

#include <filesystem>
namespace fs = std::filesystem;

#include "handlecut.h"
#include "util.h"

MAKE_FORMATTABLE(fs::path);

int main(int argc, char *argv[]) {
  SPDLOG_INFO(ANSI_BOLD "Version: {}" ANSI_RESET, VERSIONTAG);

  // Setup args
  args::ArgumentParser parser("Handle cut app");
  args::HelpFlag help(parser, "help", "Display this help message", {'h', "help"});
  args::ValueFlag<std::string> argInputFile(parser, "path", "Path to the input mesh file", {"input-file"});
  args::ValueFlag<std::string> argSeedVertices(parser, "<integer>,<integer>,...", "List of seed vertex indices for cutting", {"seed-vertices"});
  args::Flag argTunnel(parser, "tunnel", "Choose tunnel loop rather than handle loop for cutting", {"tunnel"});

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  ASSERT_WITH_LOG(argInputFile, "--input-file argument is mandatory");
  ASSERT_WITH_LOG(argSeedVertices, "--seed-vertices argument is mandatory");

  // Read input mesh
  const fs::path path_input_file = args::get(argInputFile);

  ASSERT_WITH_LOG(fs::exists(path_input_file), "File {} doesn't exist", path_input_file);
  ASSERT_WITH_LOG(path_input_file.extension() == ".obj", "Input mesh file must be in .obj format");

  MeshGeometry input;
  std::tie(input.mesh, input.geom) = readManifoldSurfaceMesh(path_input_file);

  SPDLOG_INFO("nConnectedComponents (input): {}", input.mesh->nConnectedComponents());
  SPDLOG_INFO("genus (input): {}", input.mesh->genus());

  // Parse seed vertex indices
  std::set<int> seedVertices;
  for (const std::string& s : util::tokenize(args::get(argSeedVertices), ',')) {
    seedVertices.insert(std::stoi(s));
  }

  // Run
  MeshGeometry output;
  for (int i = 0; int seedVertex : seedVertices) {
    if (i++ > 0) {
      std::swap(input.mesh, output.mesh);
      std::swap(input.geom, output.geom);
    }
    if (!handlecut(input, seedVertex, output, args::get(argTunnel)))
      SPDLOG_WARN("Error in handlecut");
  }

  SPDLOG_INFO("nConnectedComponents (output): {}", output.mesh->nConnectedComponents());
  SPDLOG_INFO("genus (output): {}", output.mesh->genus());

  // Write resulting mesh
  const fs::path path_output_file = path_input_file.parent_path() / (path_input_file.stem().string() + "-handlecut.obj");
  SPDLOG_INFO("Writing resulting mesh to {}", path_output_file);
  writeSurfaceMesh(*output.mesh, *output.geom, path_output_file);

  return 0;
}
