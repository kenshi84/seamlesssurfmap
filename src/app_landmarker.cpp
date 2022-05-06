#include <args/args.hxx>

#include <GL/glew.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl2.h>
#include <imgui_markdown.h>

#include <kt84/glfw_util.hh>
#include <kt84/geometry/CameraFree.hh>
#include <kt84/graphics/DisplayList.hh>
#include <kt84/graphics/graphics_util.hh>
#include <kt84/glut_clone/geometry.hh>
using namespace kt84::graphics_util;

#include "handlecut.h"
#include "util.h"

#include <filesystem>
namespace fs = std::filesystem;
MAKE_FORMATTABLE(fs::path);

using kt84::Camera;
using kt84::CameraFree;
using kt84::DisplayList;

//----------------+
// The model pair |
//----------------+
namespace {
struct Model : public ModelBase {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;

  AlignedBox3d bbox;
  double bboxMaxSideLength = 0.;
  CameraFree camera;
  DisplayList dispList;

  BoundaryLoopData<BoundaryLoop> correspondingBLAcrossCut;
  BoundaryLoopData<BoundaryLoop> correspondingBLAcrossModel;
  BoundaryLoopData<Vertex> handleLandmarkPerBL;

  VertexData<Vertex> correspondingLandmark;     // Only for normal landmarks
};
}
Model modelA = {"A"};
Model modelB = {"B"};

void forEachModel(std::function<void(Model& model)> f) {
  f(modelA);
  f(modelB);
};

void forEachModel2(std::function<void(Model& modelP, Model& modelQ)> f) {
  f(modelA, modelB);
  f(modelB, modelA);
};

//-----------------------------------------+
// Normal/handle landmarks from Model data |
//-----------------------------------------+
std::set<VertexPair> normalLandmarks() {
  // If uninitialized, return empty
  if (!modelA.correspondingLandmark.getMesh()) return {};

  std::set<VertexPair> res;
  for (Vertex A_v : modelA.mesh->vertices()) {
    Vertex B_v = modelA.correspondingLandmark[A_v];
    if (B_v != Vertex()) {
      ASSERT(modelB.correspondingLandmark[B_v] == A_v);
      res.insert({ A_v, B_v });
    }
  }
  return res;
}

std::set<VertexPair> handleLandmarks() {
  std::set<VertexPair> res;
  for (BoundaryLoop A_bl : modelA.mesh->boundaryLoops()) {
    BoundaryLoop B_bl = modelA.correspondingBLAcrossModel[A_bl];
    res.insert({
      modelA.handleLandmarkPerBL[A_bl],
      modelB.handleLandmarkPerBL[B_bl]
    });
  }
  return res;
}

//-----------------------------+
// Light & material parameters |
//-----------------------------+
struct LightParam {
  bool enabled = false;
  float ambient = 0.f;
  float diffuse = 0.f;
  float specular = 0.f;
  Vector4f position = {0.f, 0.f, 1.f, 0.f};
};
std::array<LightParam, 8> lightParam;

struct MaterialParam {
  float ambient = 0.2f;
  float specular = 0.5f;
  float shininess = 20.f;
} materialParam;

//------------------------+
// Other global variables |
//------------------------+
fs::path path_data_dir;
int window_width = 0;
int window_height = 0;
kt84::CameraFree* camera_active = nullptr;
ImGuiIO* io = nullptr;
bool showErrorDialog = false;
std::string errorMessage;
int genus = -1;
VertexPair vSelected;
Pair<bool> dragging = { false, false };
bool white_background = false;

//-------------------+
// Utility functions |
//-------------------+
void writeLandmarks() {
  util::write_landmarks(normalLandmarks(), handleLandmarks(), path_data_dir / "landmarks.txt");
}

BoundaryLoop boundaryVertexToBoundaryLoop(Vertex v) {
  ASSERT(v.isBoundary());
  return v.halfedge().twin().face().asBoundaryLoop();
}

Vector3d center(const AlignedBox3d& bbox) {
  return 0.5 * (bbox.min() + bbox.max());
}

//--------------------+
// Callback functions |
//--------------------+
void callback_framebuffersize(GLFWwindow* window, int width, int height) {
  window_width = width;
  window_height = height;
  modelA.camera.reshape(width, height);
  modelB.camera.reshape(width, height);
}

void mapCursorPos(GLFWwindow* window, int &x, int &y) {
  int wWidth, wHeight;
  glfwGetWindowSize(window, &wWidth, &wHeight);
  int fbWidth, fbHeight;
  glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
  x = x * fbWidth  / wWidth;
  y = y * fbHeight / wHeight;
}

void setViewportAndCameraMatrix(const Model& model);   // Forward declaration needed
void callback_mousebutton(GLFWwindow* window, int button, int action, int mods) {
  if (io->WantCaptureMouse) return;

  auto modFlag = kt84::glfw_util::parseMods(mods);
  auto actionFlag = kt84::glfw_util::parseAction(action);
  auto mouse = kt84::glfw_util::getCursorPos(window);

  mapCursorPos(window, mouse.x, mouse.y);

  if (actionFlag.press) {
    bool clickedOnLeftSide = mouse.x < window_width / 2;
    Model& model = clickedOnLeftSide ? modelA : modelB;

    // Camera control
    if (modFlag.alt) {
      camera_active = &model.camera;
      camera_active->mouse_down(mouse.x, mouse.y, modFlag.shift ? Camera::DragMode::ZOOM : modFlag.ctrl ? Camera::DragMode::PAN : Camera::DragMode::ROTATE);
      return;
    }

    float sz = read_depth(mouse.x, window_height - mouse.y);

    // Nothing happens when background clicked
    if (sz == 1.f) return;

    // Convert mouse position from 2D to 3D
    setViewportAndCameraMatrix(model);
    Vector3 mouse3D = kt84::vector_cast<3, Vector3>(unproject({(double)mouse.x, (double)(window_height - mouse.y), (double)sz}));

    // Pick vertex closest to mouse position
    kt84::MinSelector<Vertex> vClosest;
    for (Vertex v : model.mesh->vertices())
      vClosest.update(norm(model.geom->vertexPositions[v] - mouse3D), v);

    SPDLOG_INFO("Vertex closest to clicked position: {}", vClosest.value.getIndex());

    if (genus != 0) {
      // SHIFT+Click --> Specify new handle landmark
      if (modFlag.shift) {
        // Don't allow selecting boundary vertex
        if (vClosest.value.isBoundary()) return;

        model.select(vSelected) = vClosest.value;
        model.dispList.invalidate();
      }

    } else {  // genus == 0
      // SHIFT+Click --> Specify new normal landmark
      if (modFlag.shift && !modFlag.ctrl) {
        // Don't allow selecting boundary vertex
        if (vClosest.value.isBoundary()) return;

        // Don't allow selecting already landmarked vertex 
        if (model.correspondingLandmark[vClosest.value] != Vertex()) return;

        model.select(vSelected) = vClosest.value;
        model.dispList.invalidate();
        return;
      }

      // SHIFT+CTRL+Click --> Delete existing normal landmark
      if (modFlag.shift && modFlag.ctrl) {
        // Pick closest normal landmark
        vClosest = {};
        for (Vertex v : model.mesh->vertices()) {
          if (model.correspondingLandmark[v] != Vertex())
            vClosest.update(norm(model.geom->vertexPositions[v] - mouse3D), v);
        }
        if (!vClosest.count || vClosest.score > model.bboxMaxSideLength * 0.05) return;

        // Identify the selected landmark pair
        VertexPair landmark;
        if (model.name == "A") {
          landmark.A = vClosest.value;
          landmark.B = modelA.correspondingLandmark[landmark.A];
        } else {
          landmark.B = vClosest.value;
          landmark.A = modelB.correspondingLandmark[landmark.B];
        }

        // Remove it & redraw
        forEachModel([&](Model& model) {
          model.correspondingLandmark[model.select(landmark)] = Vertex();
          model.dispList.invalidate();
        });

        writeLandmarks();
        return;
      }

      // CTRL+Click --> Switch assignemnt for handle landmark
      if (!modFlag.shift && modFlag.ctrl) {
        if (!model.mesh->nBoundaryLoops()) return;    // Only relevant for meshes with handle loops

        // Pick closest handle landmark
        kt84::MinSelector<BoundaryLoop> blClosest;
        for (BoundaryLoop bl : model.mesh->boundaryLoops()) {
          if (bl < model.correspondingBLAcrossCut[bl]) {
            Vertex v = model.handleLandmarkPerBL[bl];
            blClosest.update(norm(model.geom->vertexPositions[v] - mouse3D), bl);
          }
        }
        ASSERT(blClosest.count);
        if (blClosest.score > model.bboxMaxSideLength * 0.05) return;

        // Get the other one across the cut
        BoundaryLoop blOther = model.correspondingBLAcrossCut[blClosest.value];

        // Identify the corresponding pair on the other model
        Pair<std::array<BoundaryLoop, 2>> bl;
        model.select(bl) = { blClosest.value, blOther };
        model.selectOther(bl) = {
          model.correspondingBLAcrossModel[model.select(bl)[0]],
          model.correspondingBLAcrossModel[model.select(bl)[1]]
        };

        // Swap correspondingBLAcrossModel between the two BLs across the cut
        forEachModel([&](Model& model_) {
          BoundaryLoop bl0 = model_.select(bl)[0];
          BoundaryLoop bl1 = model_.select(bl)[1];
          std::swap(
            model_.correspondingBLAcrossModel[bl0],
            model_.correspondingBLAcrossModel[bl1]
          );
          model_.dispList.invalidate();
        });

        // Sanity check
        forEachModel2([](Model& modelP, Model& modelQ) {
          for (BoundaryLoop bl : modelP.mesh->boundaryLoops())
            ASSERT(modelQ.correspondingBLAcrossModel[modelP.correspondingBLAcrossModel[bl]] == bl);
        });

        writeLandmarks();
        return;
      }

      // Click --> Start dragging for existing normal landmark
      if (!modFlag.shift && !modFlag.ctrl) {
        // Identify clicked normal landmark
        vClosest = {};
        for (Vertex v : model.mesh->vertices()) {
          // Ignore non-landmark vertices
          if (model.correspondingLandmark[v] == Vertex()) continue;

          vClosest.update(norm(model.geom->vertexPositions[v] - mouse3D), v);
        }
        if (!vClosest.count || vClosest.score > model.bboxMaxSideLength * 0.05) return;

        // Select vertex and start dragging
        model.select(dragging) = true;
        model.select(vSelected) = vClosest.value;
        model.selectOther(vSelected) = model.correspondingLandmark[vClosest.value];

        // Delete previous normal landmark
        forEachModel([](Model& model_){
          model_.correspondingLandmark[model_.select(vSelected)] = Vertex();
          model_.dispList.invalidate();
        });
      }
    }

  } else if (actionFlag.release) {
    if (camera_active) {
      camera_active->mouse_up();
      camera_active = nullptr;
      modelA.dispList.invalidate();
      modelB.dispList.invalidate();
      return;
    }

    if (dragging.A || dragging.B) {
      // Create new normal landmark
      forEachModel([](Model& model) {
        ASSERT(model.select(vSelected) != Vertex());
        model.correspondingLandmark[model.select(vSelected)] = model.selectOther(vSelected);
      });

      // Clear selection & redraw
      modelA.dispList.invalidate();
      modelB.dispList.invalidate();
      vSelected = {};
      dragging = { false, false };

      writeLandmarks();
    }
  }
}

void callback_cursorpos(GLFWwindow* window, double xpos, double ypos) {
  int mouse_x = static_cast<int>(xpos);
  int mouse_y = static_cast<int>(ypos);
  mapCursorPos(window, mouse_x, mouse_y);

  if (camera_active) {
    camera_active->mouse_move(mouse_x, mouse_y);
    return;
  }

  // Shift normal landmark by dragging
  if (dragging.A || dragging.B) {
    bool draggedOnLeftSide = mouse_x < window_width / 2;

    // If dragged mouse goes over to the other side, ignore
    if (!(dragging.A && draggedOnLeftSide || dragging.B && !draggedOnLeftSide)) return;

    Model& model = draggedOnLeftSide ? modelA : modelB;

    float sz = read_depth(mouse_x, window_height - mouse_y);

    // Nothing happens when dragging on background
    if (sz == 1.f) return;

    // Convert mouse position from 2D to 3D
    setViewportAndCameraMatrix(model);
    Vector3 mouse3D = kt84::vector_cast<3, Vector3>(unproject({(double)mouse_x, (double)(window_height - mouse_y), (double)sz}));

    // Pick vertex closest to mouse position
    kt84::MinSelector<Vertex> vClosest;
    for (Vertex v : model.mesh->vertices()) {
      // Don't allow selecting boundary vertex
      if (v.isBoundary()) continue;

      // Don't allow selecting already landmarked vertex 
      if (model.correspondingLandmark[v] != Vertex()) continue;

      vClosest.update(norm(model.geom->vertexPositions[v] - mouse3D), v);
    }

    if (model.select(vSelected) != vClosest.value)
      SPDLOG_INFO("Vertex closest to dragged position: {}", vClosest.value.getIndex());

    // Update selection & redraw
    model.select(vSelected) = vClosest.value;
    model.dispList.invalidate();
  }
}

void callback_key(GLFWwindow* window, int key, int scancode, int action, int mods) {
  auto modFlag = kt84::glfw_util::parseMods(mods);
  auto actionFlag = kt84::glfw_util::parseAction(action);
  auto mouse = kt84::glfw_util::getCursorPos(window);

  mapCursorPos(window, mouse.x, mouse.y);

  if (!actionFlag.press) return;

  // Change camera center at cursor position
  if (key == GLFW_KEY_C) {
    Model& model = mouse.x < window_width / 2 ? modelA : modelB;
    double sz = (double)read_depth(mouse.x, window_height - mouse.y);
    if (sz < 1.) {
      setViewportAndCameraMatrix(model);
      model.camera.center = unproject(Vector3d{(double)mouse.x, (double)(window_height - mouse.y), sz});
      Vector3d d = (model.camera.center - model.camera.eye).normalized();
      model.camera.up = (model.camera.up - d.dot(model.camera.up) * d).normalized();
    }
    return;
  }

  if (genus != 0) {
    // Perform handle cutting
    if (key == GLFW_KEY_ENTER) {
      if (vSelected.A != Vertex() && vSelected.B != Vertex()) {
        Pair<MeshGeometry> tempOutput;
        bool failure = false;
        VertexPair correspondingVertexAcrossCut;

        forEachModel([&](Model& model) {
          if (failure) return;

          // Annoying to have to convert this way...
          MeshGeometry tempInput = {
            std::move(model.mesh),
            std::move(model.geom)
          };

          int correspondingVertex;
          handlecut_errorcode ec;
          if (!handlecut(tempInput, model.select(vSelected).getIndex(), model.select(tempOutput), modFlag.ctrl, &correspondingVertex, &ec)) {
            failure = true;
            showErrorDialog = true;
            errorMessage = fmt::format("Handle cutting on model {} failed!\nReason:\n  {}", model.name, 
              ec == HANDLECUT_PATH_NOTFOUND ? "Loop path could not be found." :
              ec == HANDLECUT_LOOP_INVALID ? "Found loop cuts the mesh into two disconnected components." : "UNKNOWN"
            );
          }

          model.mesh = std::move(tempInput.mesh);
          model.geom = std::move(tempInput.geom);

          model.select(correspondingVertexAcrossCut) = model.select(tempOutput).mesh->vertex(correspondingVertex);
        });

        // If success, replace mesh/geom with new ones
        if (!failure) {
          Pair<std::array<BoundaryLoop, 2>> correspondingBLAcrossCut;

          // Carry over correspondingBLAcrossCut / correspondingBLAcrossModel / handleLandmarkPerBL to the new mesh
          forEachModel2([&](Model& modelP, Model& modelQ) {
            ManifoldSurfaceMesh& new_meshP = *modelP.select(tempOutput).mesh;
            ManifoldSurfaceMesh& new_meshQ = *modelQ.select(tempOutput).mesh;

            BoundaryLoopData<BoundaryLoop> new_correspondingBLAcrossCut(new_meshP);
            BoundaryLoopData<BoundaryLoop> new_correspondingBLAcrossModel(new_meshP);
            BoundaryLoopData<Vertex> new_handleLandmarkPerBL(new_meshP);

            for (BoundaryLoop bl : modelP.mesh->boundaryLoops()) {
              BoundaryLoop new_bl = new_meshP.boundaryLoop(bl.getIndex());

              new_correspondingBLAcrossCut[new_bl] = new_meshP.boundaryLoop(modelP.correspondingBLAcrossCut[bl].getIndex());
              new_correspondingBLAcrossModel[new_bl] = new_meshQ.boundaryLoop(modelP.correspondingBLAcrossModel[bl].getIndex());
              new_handleLandmarkPerBL[new_bl] = new_meshP.vertex(modelP.handleLandmarkPerBL[bl].getIndex());
            }

            modelP.correspondingBLAcrossCut = new_correspondingBLAcrossCut;
            modelP.correspondingBLAcrossModel = new_correspondingBLAcrossModel;
            modelP.handleLandmarkPerBL = new_handleLandmarkPerBL;
          });

          forEachModel([&](Model& model) {
            // Discard old mesh, replace by new
            model.mesh = std::move(model.select(tempOutput).mesh);
            model.geom = std::move(model.select(tempOutput).geom);

            model.geom->requireFaceNormals();

            // Update correspondingBLAcrossCut / handleLandmarkPerBL
            Vertex landmark0 = model.mesh->vertex(model.select(vSelected).getIndex());
            Vertex landmark1 = model.select(correspondingVertexAcrossCut);
            ASSERT(landmark0.isBoundary());
            ASSERT(landmark1.isBoundary());

            BoundaryLoop& bl0 = model.select(correspondingBLAcrossCut)[0];
            BoundaryLoop& bl1 = model.select(correspondingBLAcrossCut)[1];

            bl0 = boundaryVertexToBoundaryLoop(landmark0);
            bl1 = boundaryVertexToBoundaryLoop(landmark1);

            model.correspondingBLAcrossCut[bl0] = bl1;
            model.correspondingBLAcrossCut[bl1] = bl0;

            model.handleLandmarkPerBL[bl0] = landmark0;
            model.handleLandmarkPerBL[bl1] = landmark1;
          });

          // Clear selection & redraw
          vSelected = {};
          modelA.dispList.invalidate();
          modelB.dispList.invalidate();

          // Update correspondingBLAcrossModel
          forEachModel2([&](Model& modelP, Model& modelQ) {
            for (int i : {0, 1})
              modelP.correspondingBLAcrossModel[modelP.select(correspondingBLAcrossCut)[i]] = modelQ.select(correspondingBLAcrossCut)[i];
          });

          // When the mesh becomes genus zero, write it to *-genus0.obj
          genus = modelA.mesh->genus();
          if (genus == 0) {
            SPDLOG_INFO("The meshes are both now genus zero.");
            forEachModel([&](Model& model) {
              const fs::path path_genus0 = path_data_dir / (model.name + "-genus0.obj");
              SPDLOG_INFO("Writing to {}", path_genus0);
              writeSurfaceMesh(*model.mesh, *model.geom, path_genus0);

              // Now that genus is zero, initialize correspondingLandmark
              model.correspondingLandmark = VertexData<Vertex>(*model.mesh);
            });

            writeLandmarks();
          }
        }
      }
    }

  } else  {   // genus == 0
    if (key == GLFW_KEY_ENTER) {
      // Register the pair as landmarks
      if (vSelected.A != Vertex() && vSelected.B != Vertex()) {
        forEachModel2([&](Model& modelP, Model& modelQ) {
          modelP.correspondingLandmark[modelP.select(vSelected)] = modelQ.select(vSelected);
        });

        // Clear selection & redraw
        vSelected = {};
        modelA.dispList.invalidate();
        modelB.dispList.invalidate();

        writeLandmarks();
      }
    }
  }
}

//----------------+
// Draw functions |
//----------------+
void setViewportAndCameraMatrix(const Model& model) {
  // Draw modelA on the left half, modelB on the right half
  glViewport(&model == &modelA ? 0 : (window_width / 2), 0, window_width / 2, window_height);

  // Projection matrix
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  double zNear = model.camera.center_to_eye().norm() * 0.1;
  double zFar  = zNear * 10 + model.bboxMaxSideLength * 10;
  double aspect_ratio = window_width / 2 / (double)window_height;
  gluPerspective(40, aspect_ratio, zNear, zFar);

  // Modelview matrix
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Set light directions in the camera coordinate
  for (size_t i = 0; i < 8; ++i) {
    glLightPosition4f(lightParam[i].position, i);
  }

  gluLookAt(model.camera.get_eye(), model.camera.center, model.camera.get_up());
}

void drawVerticalGradient() {
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  // Vertical gradient
  glBegin(GL_QUADS);
  glColor3ub(0x55, 0xaa, 0xff);     // Bottom
  glVertex2d(-1, -1);
  glVertex2d(1, -1);
  glColor3d(0, 0, 0);               // Top
  glVertex2d(1, 1);
  glVertex2d(-1, 1);
  glEnd();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

void drawViewportBorder() {
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glLineWidth(5);
  glColor3d(0,0,0);
  glBegin(GL_LINE_STRIP);
  glVertex2d(-1,-1);
  glVertex2d(1,-1);
  glVertex2d(1,1);
  glVertex2d(-1,1);
  glEnd();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

void draw() {
  forEachModel([&](Model& model){
    setViewportAndCameraMatrix(model);

    if (!white_background) {
      glDepthMask(GL_FALSE);
      drawVerticalGradient();
      glDepthMask(GL_TRUE);
    }

    model.dispList.render([&](){
      // Faces with lighting
      glDisable(GL_LIGHTING);
      for (size_t i = 0; i < 8; ++i) {
        if (lightParam[i].enabled) {
          glEnable(GL_LIGHTING);
          glEnable(GL_LIGHT0 + i);
        } else {
          glDisable(GL_LIGHT0 + i);
        }
      }
      glBegin(GL_TRIANGLES);
      glColor3ub(192, 192, 192);
      for (Face f : model.mesh->faces()) {
        glNormal3d(model.geom->faceNormals[f]);
        for (Vertex v : f.adjacentVertices())
          glVertex3d(model.geom->vertexPositions[v]);
      }
      glEnd();

      auto drawSphere = [&](const Vector3& pos, double radius = 0.01) {
        glPushMatrix();
        glTranslated(pos);
        glScaled(0.7 * norm(pos - kt84::vector_cast<3, Vector3>(model.camera.get_eye())));
        kt84::glutSolidSphere(radius, 24, 12);
        glPopMatrix();
      };

      // Draw selected vertex
      if (model.select(vSelected) != Vertex()) {
        glColor3d(1, 0, 0);
        drawSphere(model.geom->vertexPositions[model.select(vSelected)], 0.012);
      }

      // Landmark vertices
      if (genus != 0) {
        // While in handle cutting, draw handleLandmarks
        for (BoundaryLoop bl : model.mesh->boundaryLoops()) {
          BoundaryLoop A_bl0 = model.name == "A" ? bl : modelB.correspondingBLAcrossModel[bl];
          BoundaryLoop A_bl1 = modelA.correspondingBLAcrossCut[A_bl0];
          if (A_bl0 > A_bl1) continue;  // Omit corresponding one across cut

          Vertex landmark = model.handleLandmarkPerBL[bl];
          Vertex A_landmark = modelA.handleLandmarkPerBL[A_bl0];

          glColor3d(util::get_random_color<int>(A_landmark.getIndex()));

          drawSphere(model.geom->vertexPositions[landmark]);
        }

      } else {
        // While in normal landmarking, draw handleLandmarks with small displacement orthogonal to boundary
        for (BoundaryLoop bl : model.mesh->boundaryLoops()) {
          BoundaryLoop A_bl = model.name == "A" ? bl : modelB.correspondingBLAcrossModel[bl];

          Vertex landmark = model.handleLandmarkPerBL[bl];
          Vertex A_landmark = modelA.handleLandmarkPerBL[A_bl];

          glColor3d(util::get_random_color<int>(A_landmark.getIndex()));

          Vector3 dirE = model.geom->halfedgeVector(landmark.halfedge());
          Vector3 dirN = model.geom->faceNormals[landmark.halfedge().face()];
          Vector3 dirO = normalize(cross(dirN, dirE));
          Vector3 offset = 0.005 * model.bboxMaxSideLength * dirO;

          drawSphere(model.geom->vertexPositions[landmark] + offset);
        }

        // Draw normalLandmarks
        for (Vertex landmark : model.mesh->vertices()) {
          if (model.correspondingLandmark[landmark] == Vertex()) continue;

          Vertex A_landmark = model.name == "A" ? landmark : modelB.correspondingLandmark[landmark];

          glColor3d(util::get_random_color<int>(A_landmark.getIndex()));

          drawSphere(model.geom->vertexPositions[landmark]);
        }
      }

      // Edges without lighting
      glDisable(GL_LIGHTING);
      glLineWidth(1.0);
      glBegin(GL_LINES);
      glColor3ub(64, 64, 64);
      for (Edge e : model.mesh->edges()) {
        for (Vertex v : e.adjacentVertices())
          glVertex3d(model.geom->vertexPositions[v]);
      }
      glEnd();
      // Boundary edges
      glLineWidth(5.0);
      glBegin(GL_LINES);
      glColor3d(0, 0, 1);
      for (Edge e : model.mesh->edges()) {
        if (e.isBoundary())
          for (Vertex v : e.adjacentVertices())
            glVertex3d(model.geom->vertexPositions[v]);
      }
      glEnd();
    });

    glDisable(GL_DEPTH_TEST);
    drawViewportBorder();
    glEnable(GL_DEPTH_TEST);
  });

  // Modal error dialog
  if (showErrorDialog) {
    ImGui::OpenPopup("Error!");
  }
  ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
  ImGui::SetNextWindowSize({1000, 1000}, ImGuiCond_Appearing);
  if (ImGui::BeginPopupModal("Error!", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
    ImGui::Markdown(errorMessage);
    ImGui::Separator();
    if (ImGui::Button("  OK  ")) {
      ImGui::CloseCurrentPopup();
      showErrorDialog = false;
      errorMessage = {};
    }
    ImGui::EndPopup();
  }

  // Camera control instructions
  ImGui::SetNextWindowPos({10, 10}, ImGuiCond_Once);
  ImGui::SetNextWindowSize({180, 230}, ImGuiCond_Once);
  // ImGui::SetNextWindowCollapsed(true, ImGuiCond_Once);
  ImGui::Begin("Camera control", nullptr, ImGuiWindowFlags_NoScrollbar);
  ImGui::Markdown(
    "  * *ALT* + *drag*\n"
    "    * Rotate\n"
    "  * *SHIFT* + *ALT* + *drag*\n"
    "    * Zoom\n"
    "  * *CTRL* + *ALT* + *drag*\n"
    "    * Pan\n"
    "  * *C* key\n"
    "    * Change center to cursor position\n"
  );
  if (ImGui::Button("Reset")) {
    forEachModel([](Model& model){
      model.camera.init(center(model.bbox) + Vector3d(0, 0, model.bboxMaxSideLength * 2.0), center(model.bbox), Vector3d::UnitY());
      model.dispList.invalidate();
    });
  }
  static char cameraParamBuf[512];
  {
    std::ostringstream oss;
    oss << " " << modelA.camera.eye.transpose();
    oss << " " << modelA.camera.center.transpose();
    oss << " " << modelA.camera.up.transpose();
    oss << " " << modelB.camera.eye.transpose();
    oss << " " << modelB.camera.center.transpose();
    oss << " " << modelB.camera.up.transpose();
    ASSERT(oss.str().size() < 512);
    strncpy(cameraParamBuf, oss.str().c_str(), oss.str().size());
    cameraParamBuf[oss.str().size()] = '\0';
  }
  if (ImGui::InputText("Param", cameraParamBuf, 512)) {
    std::vector<std::string> s;
    std::istringstream iss(cameraParamBuf);
    while (!iss.eof()) {
      s.emplace_back();
      iss >> s.back();
    }
    if (s.size() == 18) {
      std::array<double, 18> d;
      bool success = true;
      for (size_t i = 0; i < 18; ++i) {
        try {
          d[i] = std::stof(s[i]);
        } catch (...) {
          success = false;
          break;
        }
      }
      if (success) {
        size_t i = 0;
        modelA.camera.eye    << d[i], d[i + 1], d[i + 2];   i += 3;
        modelA.camera.center << d[i], d[i + 1], d[i + 2];   i += 3;
        modelA.camera.up     << d[i], d[i + 1], d[i + 2];   i += 3;
        modelB.camera.eye    << d[i], d[i + 1], d[i + 2];   i += 3;
        modelB.camera.center << d[i], d[i + 1], d[i + 2];   i += 3;
        modelB.camera.up     << d[i], d[i + 1], d[i + 2];   i += 3;
        modelA.dispList.invalidate();
        modelB.dispList.invalidate();
      }
    }
  }
  ImGui::End();

  if (genus != 0) {
    ImGui::SetNextWindowPos({10, 250}, ImGuiCond_Once);
    ImGui::SetNextWindowSize({180, 340}, ImGuiCond_Once);
    ImGui::SetNextWindowCollapsed(true, ImGuiCond_Once);
    ImGui::Begin("Handle cutting mode", nullptr, ImGuiWindowFlags_NoScrollbar);

    ImGui::Markdown(
      fmt::format("The mesh has *{}* handle(s) which you must first cut open.", genus) +
      "\n\nTo do so, *SHIFT* + *click* on each model to specify a landmark near a handle to be cut." +
      "\n\nAfter having specified the landmark on each model, press *ENTER* to perform handle cutting." +
      "\n\nUse *CTRL* + *ENTER* to choose tunnel loop rather than handle loop for cutting.");
    ImGui::End();

  } else {
    ImGui::SetNextWindowPos({10, 250}, ImGuiCond_Once);
    ImGui::SetNextWindowSize({180, 310}, ImGuiCond_Once);
    ImGui::SetNextWindowCollapsed(true, ImGuiCond_Once);
    ImGui::Begin("Landmarking mode", nullptr, ImGuiWindowFlags_NoScrollbar);

    ImGui::Markdown(
      "  * *SHIFT* + *click* to specify a landmark on each model, and *ENTER* to confirm the landmark pair.\n"
      "  * *SHIFT* + *CTRL* + *click* to delete the landmark pair.\n"
      "  * *Drag* to shift landmark to nearby vertex.\n"
      "  * *CTRL* + *click* on handle landmarks to switch between the two possible pair assignments.\n");
    ImGui::End();
  }

}

//------+
// Main |
//------+
int main(int argc, char *argv[]) {
  SPDLOG_INFO(ANSI_BOLD "Version: {}" ANSI_RESET, VERSIONTAG);

  // Setup args
  args::ArgumentParser parser("Landmarker");
  args::HelpFlag help(parser, "help", "Display this help message", {'h', "help"});
  args::ValueFlag<std::string> argDataDir(parser, "path", "Path to the data directory to read from / write to (default: current directory)", {"data-dir"});
  args::Flag argWhiteBackground(parser, "white-background", "Use white background instead of vertical color gradient", {"white-background"});
  args::Flag argWriteNormalized(parser, "write-normalized", "Write area-normalized mesh as *-normalized.obj", {"write-normalized"});

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
  path_data_dir = argDataDir ? fs::path(args::get(argDataDir)) : fs::current_path();
  ensure_existence(path_data_dir);

  white_background = args::get(argWhiteBackground);

  SPDLOG_INFO("Data directory: {}", path_data_dir);

  // Read orig mesh & topology sanity check
  forEachModel([&](Model& model) {
    const fs::path path_orig = path_data_dir / (model.name + "-orig.obj");
    ensure_existence(path_orig);

    std::tie(model.mesh, model.geom) = readManifoldSurfaceMesh(path_orig);

    // Compute bouding box
    for (Vertex v : model.mesh->vertices())
      model.bbox.extend(kt84::vector_cast<3, Vector3d>(model.geom->vertexPositions[v]));

    model.bboxMaxSideLength = model.bbox.diagonal().maxCoeff();

    model.camera.init(center(model.bbox) + Vector3d(0, 0, model.bboxMaxSideLength * 2.0), center(model.bbox), Vector3d::UnitY());

    model.geom->requireFaceNormals();

    ASSERT_WITH_LOG(model.mesh->nConnectedComponents() == 1, "{}-orig has multiple connected components ({})", model.name, model.mesh->nConnectedComponents());
    ASSERT_WITH_LOG(model.mesh->nBoundaryLoops() == 0, "{}-orig is not closed", model.name);

    SPDLOG_INFO("{}'s' nVertices is {}, nFaces is {}", model.name, model.mesh->nVertices(), model.mesh->nFaces());

    // Write area-normalized mesh as *-normalized
    if (args::get(argWriteNormalized)) {
      std::unique_ptr<VertexPositionGeometry> geom_normalized = model.geom->copy();

      geom_normalized->requireFaceAreas();

      double surfaceArea = 0;
      for (Face f : model.mesh->faces()) {
        surfaceArea += geom_normalized->faceAreas[f];
      }

      for (Vertex v : model.mesh->vertices()) {
        geom_normalized->inputVertexPositions[v] /= std::sqrt(surfaceArea);
      }

      SPDLOG_INFO("Writing area-normalized mesh to {}-normalized.obj", model.name);
      writeSurfaceMesh(*model.mesh, *geom_normalized, path_data_dir / (model.name + "-normalized.obj"));
    }
  });
  genus = modelA.mesh->genus();
  ASSERT_WITH_LOG(genus == modelB.mesh->genus(), "Genus of orig pair doesn't match ({} vs {})", genus, modelB.mesh->genus());
  SPDLOG_INFO("Surface genus is {}", genus);

  // Read landmarks if exists
  std::set<IntPair> landmarks_i;
  if (fs::exists(path_data_dir / "landmarks.txt")) {
    SPDLOG_INFO(ANSI_BOLD "Found landmark.txt, so we try to read it." ANSI_RESET);
    SPDLOG_INFO(ANSI_BOLD "In case of failure, remove the file and re-run the program." ANSI_RESET);

    std::ifstream ifs(path_data_dir / "landmarks.txt");

    std::string landmarks_str;
    std::getline(ifs, landmarks_str);

    landmarks_i = util::parse_landmarks(landmarks_str);
  }

  // If genus zero, initialize correspondingLandmark
  if (genus == 0) {
    forEachModel([](Model& model) {
      model.correspondingLandmark = VertexData<Vertex>(*model.mesh);
    });

    for (const IntPair& landmark_i : landmarks_i) {
      VertexPair landmark;
      forEachModel([&](Model& model){
        // Index bound check
        ASSERT_WITH_LOG(model.select(landmark_i) < model.mesh->nVertices(), "{}'s landmark {} is no less than nVertices ({})", model.name, model.select(landmark_i), model.mesh->nVertices());
        model.select(landmark) = model.mesh->vertex(model.select(landmark_i));
      });

      modelA.correspondingLandmark[landmark.A] = landmark.B;
      modelB.correspondingLandmark[landmark.B] = landmark.A;
    }

  } else do {
    // For nonzero genus, check if *-orig.obj and landmarks.txt exist and are proper
    if (!fs::exists(path_data_dir / "landmarks.txt")) break;
    if (!fs::exists(path_data_dir / "A-genus0.obj")) break;
    if (!fs::exists(path_data_dir / "B-genus0.obj")) break;

    SPDLOG_INFO(ANSI_BOLD "Found A-genus0.obj and B-genus0.obj along with landmarks.txt, so we try to read them and perform consistency check." ANSI_RESET);
    SPDLOG_INFO(ANSI_BOLD "In case of failure, remove these files and re-run the program." ANSI_RESET);

    // Read *-genus0.obj
    forEachModel([&](Model& model){
      MeshGeometry genus0;
      std::tie(genus0.mesh, genus0.geom) = readManifoldSurfaceMesh(path_data_dir / (model.name + "-genus0.obj"));

      ASSERT_WITH_LOG(util::compare_orig_vs_cut(*model.mesh, *model.geom, *genus0.mesh, *genus0.geom), "{}-orig and {}-genus0 don't coincide exactly", model.name, model.name);
      ASSERT_WITH_LOG(genus0.mesh->genus() == 0, "Genus of {}-genus0.obj is actually not zero", model.name);

      // Replace model.mesh/geom with genus0
      model.mesh = std::move(genus0.mesh);
      model.geom = std::move(genus0.geom);

      model.geom->requireFaceNormals();

      // Initialize correspondingBLAcrossCut & correspondingBLAcrossModel & handleLandmarkPerBL & correspondingLandmark
      model.correspondingBLAcrossCut = BoundaryLoopData<BoundaryLoop>(*model.mesh);
      model.correspondingBLAcrossModel = BoundaryLoopData<BoundaryLoop>(*model.mesh);
      model.handleLandmarkPerBL = BoundaryLoopData<Vertex>(*model.mesh);

      model.correspondingLandmark = VertexData<Vertex>(*model.mesh);
    });

    genus = 0;

    // Fill correspondingBLAcrossCut & correspondingBLAcrossModel & handleLandmarkPerBL & correspondingLandmark
    while(!landmarks_i.empty()) {
      // Pop the front element
      auto landmark0_i = landmarks_i.begin();
      landmarks_i.erase(landmark0_i);

      // Convert it from int-based to element-based
      VertexPair landmark0;
      forEachModel([&](Model& model){
        // Index bound check
        ASSERT_WITH_LOG(model.select(*landmark0_i) < model.mesh->nVertices(), "{}'s landmark {} is no less than nVertices ({})", model.name, model.select(*landmark0_i), model.mesh->nVertices());
        model.select(landmark0) = model.mesh->vertex(model.select(*landmark0_i));
      });

      // Check consistency
      ASSERT_WITH_LOG(landmark0.A.isBoundary() == landmark0.B.isBoundary(), "A's landmark {} and B's landmark {} are not simultaneously boundary/non-boundary", landmark0.A.getIndex(), landmark0.B.getIndex());

      if (landmark0.A.isBoundary()) {
        // If on boundary, there must exist corresponding boundary landmark that coincides exactly
        kt84::MinSelector<IntPair> landmark1_i;
        for (const IntPair& landmark_i : landmarks_i) {
          double dist = 0;
          forEachModel([&](Model& model){
            // Index bound check
            ASSERT_WITH_LOG(model.select(landmark_i) < model.mesh->nVertices(), "{}'s landmark {} is no less than nVertices ({})", model.name, model.select(landmark_i), model.mesh->nVertices());

            Vector3 p0 = model.geom->vertexPositions[model.select(landmark0)];
            Vector3 p1 = model.geom->vertexPositions[model.mesh->vertex(model.select(landmark_i))];

            dist += norm(p1 - p0);
          });
          landmark1_i.update(dist, landmark_i);
        }
        ASSERT_WITH_LOG(landmark1_i.score == 0.0, "Handle landmark {} must have another exactly coinciding handle landmark, but was not found", landmark0);

        // Remove it from the set
        landmarks_i.erase(landmark1_i.value);

        // Convert it from int-based to element-based
        VertexPair landmark1 = {
            modelA.mesh->vertex(landmark1_i.value.A),
            modelB.mesh->vertex(landmark1_i.value.B),
        };

        // Check consistency
        forEachModel([&](Model& model){
          ASSERT_WITH_LOG(model.select(landmark1).isBoundary(), "{}'s landmark {} is supposed to be on boundary, but not", model.name, model.select(landmark1));
        });

        // Now with landmark0 & landmark1 at hand, register them with correspondingBLAcrossCut & correspondingBLAcrossModel & handleLandmarkPerBL
        // while checking consistency
        forEachModel([&](Model& model){
          BoundaryLoop bl0 = boundaryVertexToBoundaryLoop(model.select(landmark0));
          BoundaryLoop bl1 = boundaryVertexToBoundaryLoop(model.select(landmark1));

          // Uniqueness check
          ASSERT(model.correspondingBLAcrossCut[bl0] == BoundaryLoop());
          ASSERT(model.correspondingBLAcrossCut[bl1] == BoundaryLoop());
          ASSERT(model.correspondingBLAcrossModel[bl0] == BoundaryLoop());
          ASSERT(model.correspondingBLAcrossModel[bl1] == BoundaryLoop());
          ASSERT(model.handleLandmarkPerBL[bl0] == Vertex());
          ASSERT(model.handleLandmarkPerBL[bl1] == Vertex());

          // Register
          model.correspondingBLAcrossCut[bl0] = bl1;
          model.correspondingBLAcrossCut[bl1] = bl0;
          model.correspondingBLAcrossModel[bl0] = boundaryVertexToBoundaryLoop(model.selectOther(landmark0));
          model.correspondingBLAcrossModel[bl1] = boundaryVertexToBoundaryLoop(model.selectOther(landmark1));
          model.handleLandmarkPerBL[bl0] = model.select(landmark0);
          model.handleLandmarkPerBL[bl1] = model.select(landmark1);
        });

      } else {
        // Normal landmark --> register with correspondingLandmark
        modelA.correspondingLandmark[landmark0.A] = landmark0.B;
        modelB.correspondingLandmark[landmark0.B] = landmark0.A;
      }
    }

  } while (false);

  kt84::glfw_util::InitConfig initConfig;
  initConfig.window.title = "Landmarker";
  initConfig.callback.framebuffersize = callback_framebuffersize;
  initConfig.callback.mousebutton = callback_mousebutton;
  initConfig.callback.cursorpos = callback_cursorpos;
  initConfig.callback.key = callback_key;

  kt84::glfw_util::LoopControl loopControl = kt84::glfw_util::init(initConfig);
  loopControl.idle_func[0] = []() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Set light & material parameters
    for (size_t i = 0; i < 8; ++i) {
      glLightAmbient3f (Vector3::constant(lightParam[i].ambient ), i);
      glLightDiffuse3f (Vector3::constant(lightParam[i].diffuse ), i);
      glLightSpecular3f(Vector3::constant(lightParam[i].specular), i);
    }
    glMaterialAmbient3f(Vector3::constant(materialParam.ambient));
    glMaterialSpecular3f(Vector3::constant(materialParam.specular));
    glMaterialShininess1f(materialParam.shininess);

    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    if (io->WantCaptureMouse && io->MouseReleased[0]) {
      modelA.dispList.invalidate();
      modelB.dispList.invalidate();
    }
  };
  loopControl.idle_func[1] = draw;
  loopControl.idle_func[2] = []() {
    ImGui::Render();
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
  };


  if (glewInit() != GLEW_OK) {
    SPDLOG_ERROR("Failed to initialize OpenGL context");
    return 1;
  }

  // OpenGL initial settings
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);             // Automatic normal normalization after scaling
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1, 1);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glClearColor(1, 1, 1, 1);

  // Light parameters
  lightParam[0].enabled = true;
  lightParam[0].ambient = 0.2;
  lightParam[0].diffuse = 1.0;
  lightParam[0].specular = 0.0;
  lightParam[0].position = Vector4f{ -0.40464258, 0.59512496, 0.6943284, 0.0 };
  lightParam[1].enabled = true;
  lightParam[1].ambient = 0.0;
  lightParam[1].diffuse = 0.5;
  lightParam[1].specular = 0.0;
  lightParam[1].position = Vector4f{ 0.8579609, -0.21099381, 0.46838602, 0.0 };
  lightParam[2].enabled = true;
  lightParam[2].ambient = 0.0;
  lightParam[2].diffuse = 0.25;
  lightParam[2].specular = 0.0;
  lightParam[2].position = Vector4f{ -0.5633799, -0.8006998, 0.2036816, 0.0 };

  // Material parameters
  materialParam.ambient = 0.1;
  materialParam.specular = 0.0;
  materialParam.shininess = 0.0;

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  io = &ImGui::GetIO();
  io->IniFilename = nullptr;  // Suppress generating imgui.ini
  ImGui_ImplGlfw_InitForOpenGL(loopControl.window, true);
  ImGui_ImplOpenGL2_Init();

  kt84::glfw_util::start_loop(loopControl);

  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  return 0;
}
