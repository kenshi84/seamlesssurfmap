#include <iostream>

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
using namespace kt84::graphics_util;

#include "common.h"

using kt84::Camera;
using kt84::CameraFree;
using kt84::DisplayList;

//-----------------+
// Mesh & geometry |
//-----------------+
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geom;

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
std::string filename;
AlignedBox3d bbox;
double bboxMaxSideLength = 0.;
CameraFree camera;
DisplayList dispList;
int window_width = 0;
int window_height = 0;
ImGuiIO* io = nullptr;

//-------------------+
// Utility functions |
//-------------------+
Vector3d center(const AlignedBox3d& bbox) {
  return 0.5 * (bbox.min() + bbox.max());
}

//--------------------+
// Callback functions |
//--------------------+
void callback_framebuffersize(GLFWwindow* window, int width, int height) {
  window_width = width;
  window_height = height;
  camera.reshape(width, height);
}

void mapCursorPos(GLFWwindow* window, int &x, int &y) {
  int wWidth, wHeight;
  glfwGetWindowSize(window, &wWidth, &wHeight);
  int fbWidth, fbHeight;
  glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
  x = x * fbWidth  / wWidth;
  y = y * fbHeight / wHeight;
}

void callback_mousebutton(GLFWwindow* window, int button, int action, int mods) {
  if (io->WantCaptureMouse) return;

  auto modFlag = kt84::glfw_util::parseMods(mods);
  auto actionFlag = kt84::glfw_util::parseAction(action);
  auto mouse = kt84::glfw_util::getCursorPos(window);

  mapCursorPos(window, mouse.x, mouse.y);

  // Camera control
  if (actionFlag.press) {
    camera.mouse_down(mouse.x, mouse.y, modFlag.shift ? Camera::DragMode::ZOOM : modFlag.ctrl ? Camera::DragMode::PAN : Camera::DragMode::ROTATE);

  } else if (actionFlag.release) {
    camera.mouse_up();
  }
}

void callback_cursorpos(GLFWwindow* window, double xpos, double ypos) {
  int mouse_x = static_cast<int>(xpos);
  int mouse_y = static_cast<int>(ypos);
  mapCursorPos(window, mouse_x, mouse_y);

  if (camera.drag_mode != Camera::DragMode::NONE) {
    camera.mouse_move(mouse_x, mouse_y);
  }
}

void callback_key(GLFWwindow* window, int key, int scancode, int action, int mods) {
  auto modFlag = kt84::glfw_util::parseMods(mods);
  auto actionFlag = kt84::glfw_util::parseAction(action);
  auto mouse = kt84::glfw_util::getCursorPos(window);

  mapCursorPos(window, mouse.x, mouse.y);

  if (!actionFlag.press) return;

  // Project & write
  if (key == GLFW_KEY_ENTER) {
    // Get UV per vertex
    VertexData<Vector2d> uvPerVertex(*mesh);
    AlignedBox2d uvBbox;
    for (Vertex v : mesh->vertices()) {
       uvPerVertex[v] = project(kt84::vector_cast<3, Vector3d>(geom->vertexPositions[v])).head(2);
       uvBbox.extend(uvPerVertex[v]);
    }

    // Normalize UV to [0, 1] x [0, 1]
    double uvBboxMaxSideLength = uvBbox.diagonal().maxCoeff();
    for (Vertex v : mesh->vertices()) {
      uvPerVertex[v] -= uvBbox.min();
      uvPerVertex[v] /= uvBboxMaxSideLength;
    }

    // Write texutred mesh as obj
    CornerData<Vector2> texCoords(*mesh);
    for (Corner c : mesh->corners())
      texCoords[c] = kt84::vector_cast<2, Vector2>(uvPerVertex[c.vertex()]);

    writeSurfaceMesh(*mesh, *geom, texCoords, filename + "-uv.obj");
    SPDLOG_INFO("Mesh with UV written to {}", filename + "-uv.obj");
  }
}

//----------------+
// Draw functions |
//----------------+
void setViewportAndCameraMatrix() {
  glViewport(0, 0, window_width, window_height);

  // Projection matrix
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  double zNear = camera.center_to_eye().norm() * 0.1;
  double zFar  = zNear * 10 + bboxMaxSideLength * 10;
  double aspect_ratio = window_width / (double)window_height;

  double ortho_y = zNear * 2.5;
  double ortho_x  = ortho_y * aspect_ratio;
  glOrtho(-ortho_x, ortho_x, -ortho_y, ortho_y, zNear, zFar);

  // Modelview matrix
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Set light directions in the camera coordinate
  for (size_t i = 0; i < 8; ++i) {
    glLightPosition4f(lightParam[i].position, i);
  }

  gluLookAt(camera.get_eye(), camera.center, camera.get_up());
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

void draw() {
  setViewportAndCameraMatrix();

  glDepthMask(GL_FALSE);
  drawVerticalGradient();
  glDepthMask(GL_TRUE);

  dispList.render([&](){
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
    for (Face f : mesh->faces()) {
      glNormal3d(geom->faceNormals[f]);
      for (Vertex v : f.adjacentVertices())
        glVertex3d(geom->vertexPositions[v]);
    }
    glEnd();

    // Edges without lighting
    glDisable(GL_LIGHTING);
    glLineWidth(1.0);
    glBegin(GL_LINES);
    glColor3ub(64, 64, 64);
    for (Edge e : mesh->edges()) {
      for (Vertex v : e.adjacentVertices())
        glVertex3d(geom->vertexPositions[v]);
    }
    glEnd();
  });

  // Camera control instructions
  ImGui::SetNextWindowPos({10, 10}, ImGuiCond_Once);
  ImGui::SetNextWindowSize({180, 300}, ImGuiCond_Once);
  ImGui::SetNextWindowCollapsed(true, ImGuiCond_Once);
  ImGui::Begin("Usage", nullptr, ImGuiWindowFlags_NoScrollbar);
  ImGui::Markdown(
    "# Camera control\n"
    "  * *drag*\n"
    "    * Rotate\n"
    "  * *SHIFT* + *drag*\n"
    "    * Zoom\n"
    "  * *CTRL* + *drag*\n"
    "    * Pan\n"
    "# Instruction\n"
    "  * Press *ENTER* to write mesh with UVs\n"
  );
  ImGui::End();
}

//------+
// Main |
//------+
int main(int argc, char *argv[]) {
  // Setup args
  args::ArgumentParser parser("UV projector");
  args::HelpFlag help(parser, "help", "Display this help message", {'h', "help"});
  args::Positional<std::string> argMesh(parser, "mesh", "The mesh file");

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

  if (!argMesh) {
    SPDLOG_ERROR("Please specify the mesh file.");
    return 1;
  }

  filename = args::get(argMesh);

  std::tie(mesh, geom) = readManifoldSurfaceMesh(filename);

  // Compute bouding box
  for (Vertex v : mesh->vertices())
    bbox.extend(kt84::vector_cast<3, Vector3d>(geom->vertexPositions[v]));

  bboxMaxSideLength = bbox.diagonal().maxCoeff();

  camera.init(center(bbox) + Vector3d(0, 0, bboxMaxSideLength * 2.0), center(bbox), Vector3d::UnitY());

  geom->requireFaceNormals();

  SPDLOG_INFO("{} vertices, {} faces, genus {}", mesh->nVertices(), mesh->nFaces(), mesh->genus());

  kt84::glfw_util::InitConfig initConfig;
  initConfig.window.title = "UV Projector";
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
