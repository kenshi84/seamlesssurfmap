#include <args/args.hxx>

#include <filesystem>

namespace fs = std::filesystem;

#include "util.h"
#include "cocut.h"
#include "gmap.h"
#include "lift.h"
#include "refine.h"

MAKE_FORMATTABLE(fs::path);

namespace {

struct Model : public ModelBase {
  MeshGeometry orig;                // Original closed mesh without any modification
  MeshGeometry genus0;              // For high genus case, genus zero mesh obtained after handle cutting
  MeshGeometry disk;                // Disk topology mesh after cutting genus0 along cut graph connecting landmarks

  std::unique_ptr<VertexPositionGeometry> gmap_geom;      // Flattened version of disk obtained from G-mapping

  VertexData<SurfacePoint> vertexImage;
};

}

int main(int argc, char *argv[]) {
  SPDLOG_INFO(ANSI_BOLD "Version: {}" ANSI_RESET, VERSIONTAG);

  // Declare the model pair
  Model modelA = {"A"};
  Model modelB = {"B"};

  auto forEachModel = [&](std::function<void(Model& model)> f) {
    f(modelA);
    f(modelB);
  };

  auto forEachModel2 = [&](std::function<void(Model& modelP, Model& modelQ)> f) {
    f(modelA, modelB);
    f(modelB, modelA);
  };

  // Setup args
  args::ArgumentParser parser("Seamless Surface Mappings");
  args::HelpFlag help(parser, "help", "Display this help message", {'h', "help"});
  args::ValueFlag<std::string> argDataDir(parser, "path", "Path to the data directory to read from / write to (default: current directory)", {"data-dir"});
  args::ValueFlag<int> argNumIter(parser, "integer", "Number of iterations for G-mapping optimization. Negative value means to iterate indefinitely (default: -1).", {"num-iter"});
  args::Flag argDebugOutput(parser, "debugOutput", "Enable outputting debug files", {"debug-output"});
  args::Flag argDontRefine(parser, "dontRefine", "Do not perform triangle refinement around landmarks and between close-by seams (increases chance of failure)", {"dont-refine"});

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

  auto ensure_existence = [](const fs::path& p) {
    ASSERT_WITH_LOG(fs::exists(p), "Path {} doesn't exist", p);
  };

  // Set data directory
  const fs::path path_data_dir = argDataDir ? fs::path(args::get(argDataDir)) : fs::current_path();
  ensure_existence(path_data_dir);

  SPDLOG_INFO("Data directory: {}", path_data_dir);

  // Read orig mesh & topology sanity check
  forEachModel([&](Model& model) {
    const fs::path path_orig = path_data_dir / (model.name + "-orig.obj");
    ensure_existence(path_orig);

    std::tie(model.orig.mesh, model.orig.geom, model.orig.uv) = readParameterizedManifoldSurfaceMesh(path_orig);

    ASSERT_WITH_LOG(model.orig.mesh->nConnectedComponents() == 1, "{}-orig has multiple connected components ({})", model.name, model.orig.mesh->nConnectedComponents());
    ASSERT_WITH_LOG(model.orig.mesh->nBoundaryLoops() == 0, "{}-orig is not closed", model.name);
  });

  const int nHandles = modelA.orig.mesh->genus();

  ASSERT_WITH_LOG(modelB.orig.mesh->genus() == nHandles, "Genus of orig pair doesn't match ({} vs {})", nHandles, modelB.orig.mesh->genus());

  // For nonzero genus, load genus zero meshes produced by handle cutting
  if (nHandles > 0) {
    SPDLOG_INFO("Nonzero genus case, loading genus0 meshes");

    // Read genus0 mesh & topology sanity check
    forEachModel([&](Model& model) {
      const fs::path path_genus0 = path_data_dir / (model.name + "-genus0.obj");
      ensure_existence(path_genus0);

      std::tie(model.genus0.mesh, model.genus0.geom) = readManifoldSurfaceMesh(path_genus0);

      ASSERT_WITH_LOG(model.genus0.mesh->nConnectedComponents() == 1, "{}-genus0 has multiple connected components ({})", model.name, model.genus0.mesh->nConnectedComponents());
      ASSERT_WITH_LOG(model.genus0.mesh->genus() == 0, "{}-genus0's genus is not zero ({})", model.name, model.genus0.mesh->genus());
      ASSERT_WITH_LOG(model.genus0.mesh->nBoundaryLoops() == 2 * nHandles, "{}-genus0 has wrong number of boundary loops ({})", model.name, model.genus0.mesh->nBoundaryLoops());

      ASSERT_WITH_LOG(util::compare_orig_vs_cut(model.orig, model.genus0), "{}-orig and {}-genus0 don't coincide exactly", model.name, model.name);
    });
  }

  // Read landmarks
  std::set<IntPair> landmarks;
  {
    const fs::path path_landmarks = path_data_dir / "landmarks.txt";
    ensure_existence(path_landmarks);

    std::ifstream ifs(path_landmarks);

    std::string landmarks_str;
    std::getline(ifs, landmarks_str);

    landmarks = util::parse_landmarks(landmarks_str);
  }



  SPDLOG_INFO(ANSI_BOLD "Step 1: Co-cutting" ANSI_RESET);

  std::vector<std::set<IntPair>> normalLandmarkGroups;
  std::vector<std::set<IntPair>> handleLandmarkGroups;

  cocut(
    nHandles == 0 ? modelA.orig : modelA.genus0,
    nHandles == 0 ? modelB.orig : modelB.genus0,
    landmarks,
    modelA.disk,
    modelB.disk,
    normalLandmarkGroups,
    handleLandmarkGroups,
    args::get(argDebugOutput));

  forEachModel([&](Model& model) {
    SPDLOG_INFO("Writing resulting mesh to {}-disk.obj", model.name);
    writeSurfaceMesh(*model.disk.mesh, *model.disk.geom, path_data_dir / (model.name + "-disk.obj"));
  });



  if (!args::get(argDontRefine)) {
    SPDLOG_INFO(ANSI_BOLD "Step 1.5: Refinement around landmarks" ANSI_RESET);

    forEachModel([&](Model& model){
      // Aggregate all landmarks
      std::set<int> landmarks;
      for (const std::vector<std::set<IntPair>>& landmarkGroups : { normalLandmarkGroups, handleLandmarkGroups }) {
        for (const std::set<IntPair>& landmarkGroup : landmarkGroups)
          for (const IntPair& landmark : landmarkGroup)
            landmarks.insert(model.select(landmark));
      }

      refine(landmarks, model.orig, model.genus0, model.disk);

      // Write results
      if (model.orig.uv) {
        writeSurfaceMesh(*model.orig.mesh, *model.orig.geom, *model.orig.uv, path_data_dir / (model.name + "-orig-refined.obj"));
      } else {
        writeSurfaceMesh(*model.orig.mesh, *model.orig.geom, path_data_dir / (model.name + "-orig-refined.obj"));
      }
      if (nHandles > 0) {
        writeSurfaceMesh(*model.genus0.mesh, *model.genus0.geom, path_data_dir / (model.name + "-genus0-refined.obj"));
      }
      writeSurfaceMesh(*model.disk.mesh, *model.disk.geom, path_data_dir / (model.name + "-disk-refined.obj"));
    });
  }



  SPDLOG_INFO(ANSI_BOLD "Step 2: G-mapping optimization" ANSI_RESET);

  int num_iter = argNumIter ? args::get(argNumIter) : -1;
  if (num_iter < 0) num_iter = std::numeric_limits<int>::max();

  gmap(
    *modelA.disk.mesh,
    *modelB.disk.mesh,
    *modelA.disk.geom,
    *modelB.disk.geom,
    normalLandmarkGroups,
    handleLandmarkGroups,
    num_iter,
    modelA.gmap_geom,
    modelB.gmap_geom,
    args::get(argDebugOutput));

  forEachModel([&](Model& model) {
    SPDLOG_INFO("Writing resulting mesh to {}-gmap.obj", model.name);
    writeSurfaceMesh(*model.disk.mesh, *model.gmap_geom, path_data_dir / (model.name + "-gmap.obj"));
  });



  SPDLOG_INFO(ANSI_BOLD "Step 3: Path lifting" ANSI_RESET);

  MeshGeometry A_gmap = { std::move(modelA.disk.mesh), std::move(modelA.gmap_geom) };
  MeshGeometry B_gmap = { std::move(modelB.disk.mesh), std::move(modelB.gmap_geom) };

  lift(
    modelA.orig,
    modelB.orig,
    A_gmap,
    B_gmap,
    normalLandmarkGroups,
    handleLandmarkGroups,
    modelA.vertexImage,
    modelB.vertexImage);

  forEachModel2([&](Model& modelP, Model& modelQ) {
    SPDLOG_INFO("Writing vertex image to {}to{}.txt", modelP.name, modelQ.name);

    std::ofstream ofs(path_data_dir / (modelP.name + "to" + modelQ.name + ".txt"));
    for (Vertex vP : modelP.orig.mesh->vertices()) {
      SurfacePoint spQ = modelP.vertexImage[vP];
      for (int i = 0; Vertex vQ : spQ.face.adjacentVertices()) {
        if (spQ.faceCoords[i] > 0.0) {
          ofs << vP.getIndex() << " " << vQ.getIndex() << " " << spQ.faceCoords[i] << "\n";
        }
        ++i;
      }
    }

    SPDLOG_INFO("Writing deformed mesh to {}on{}.obj", modelP.name, modelQ.name);

    VertexPositionGeometry P_geom(*modelP.orig.mesh);
    for (Vertex vP : modelP.orig.mesh->vertices()) {
      P_geom.inputVertexPositions[vP] = modelP.vertexImage[vP].interpolate(modelQ.orig.geom->vertexPositions);
    }
    writeSurfaceMesh(*modelP.orig.mesh, P_geom, path_data_dir / (modelP.name + "on" + modelQ.name + ".obj"));
  });

  forEachModel2([&](Model& modelP, Model& modelQ) {
    if (modelP.orig.uv) {
      SPDLOG_INFO("Writing {}-with-{}-uv.obj", modelQ.name, modelP.name);

      // Construct per-vertex UV on P
      VertexData<Vector2> P_uv(*modelP.orig.mesh);
      for (Vertex vP : modelP.orig.mesh->vertices())
        P_uv[vP] = (*modelP.orig.uv)[vP.corner()];

      // Transfer P_uv to Q
      CornerData<Vector2> Q_uv(*modelQ.orig.mesh);
      for (Vertex vQ : modelQ.orig.mesh->vertices()) {
        for (Corner cQ : vQ.adjacentCorners())
          Q_uv[cQ] = modelQ.vertexImage[vQ].interpolate(P_uv);
      }
      writeSurfaceMesh(*modelQ.orig.mesh, *modelQ.orig.geom, Q_uv, path_data_dir / (modelQ.name + "-with-" + modelP.name + "-uv.obj"));
    }
  });

  return 0;
}
