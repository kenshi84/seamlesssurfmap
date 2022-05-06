cmake_minimum_required(VERSION 3.1)

project(seamlesssurfmap)

# Print the build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

#--------------------+
# External libraries |
#--------------------+

include_directories(${CMAKE_CURRENT_LIST_DIR}/external)

# MOSEK
find_package(MOSEK REQUIRED)
include_directories(${MOSEK_INCLUDE_DIR})
function (fix_mosek_lib_path app)
  # Convince cmake to fix the linker path for mosek
  # http://www.alecjacobson.com/weblog/?p=4878
  add_custom_command(TARGET ${app} POST_BUILD COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change `otool -L ${app} | sed -n -e \""s/.*\\(libmosek.*dylib\\).*/\\1/p"\"` ${MOSEK_C_LIBRARY} ${app})
  add_custom_command(TARGET ${app} POST_BUILD COMMAND ${CMAKE_INSTALL_NAME_TOOL} -change `otool -L ${app} | sed -n -e \""s/.*\\(libfusion.*dylib\\).*/\\1/p"\"` ${MOSEK_CXX_LIBRARY} ${app})
endfunction()

# geometry-central
set(GC_EIGEN_LOCATION ${CMAKE_CURRENT_LIST_DIR}/external/eigen CACHE PATH "Eigen path")
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/external/geometry-central)
include_directories(${CMAKE_CURRENT_LIST_DIR}/external/geometry-central/include)

# spdlog
include_directories(${CMAKE_CURRENT_LIST_DIR}/external/spdlog/include)

# libigl
include_directories(${CMAKE_CURRENT_LIST_DIR}/external/libigl/include)
add_definitions(-DLIBIGL_TUTORIAL_DATA=\"${CMAKE_CURRENT_LIST_DIR}/external/libigl-tutorial-data\")

# UltiMaille
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/external/ultimaille)
include_directories(${CMAKE_CURRENT_LIST_DIR}/external/ultimaille)
include_directories(${CMAKE_CURRENT_LIST_DIR}/external/ultimaille/ext)

# Eigen
include_directories(${CMAKE_CURRENT_LIST_DIR}/external/eigen)

# GLEW
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/external/glew-2.2.0/build/cmake)
include_directories(${CMAKE_CURRENT_LIST_DIR}/external/glew-2.2.0/include)

## GLFW
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/external/glfw)

# ImGUI (CMakeLists.txt not provided)
include_directories(${CMAKE_CURRENT_LIST_DIR}/external/imgui)
include_directories(${CMAKE_CURRENT_LIST_DIR}/external/imgui/backends)
add_library(imgui
  ${CMAKE_CURRENT_LIST_DIR}/external/imgui/imgui.cpp
  # ${CMAKE_CURRENT_LIST_DIR}/external/imgui/imgui_demo.cpp       # Useful for testing via ImGui::ShowDemoWindow
  ${CMAKE_CURRENT_LIST_DIR}/external/imgui/imgui_draw.cpp
  ${CMAKE_CURRENT_LIST_DIR}/external/imgui/imgui_tables.cpp
  ${CMAKE_CURRENT_LIST_DIR}/external/imgui/imgui_widgets.cpp
  ${CMAKE_CURRENT_LIST_DIR}/external/imgui/backends/imgui_impl_glfw.cpp
  ${CMAKE_CURRENT_LIST_DIR}/external/imgui/backends/imgui_impl_opengl2.cpp)
target_link_libraries(imgui glfw)

#---------+
# My code |
#---------+
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

# hard-coded data path
add_definitions(-DDATA_PATH=\"${CMAKE_CURRENT_LIST_DIR}/data\")

# Version
find_package(Git QUIET)
set(VERSIONTAG "unknown")
if (GIT_FOUND OR Git_FOUND)
  execute_process(COMMAND "${GIT_EXECUTABLE}" rev-parse --short HEAD RESULT_VARIABLE RET OUTPUT_VARIABLE COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (NOT RET)
    set(VERSIONTAG "${COMMIT}")
  endif()
endif()
message(STATUS "Version tag: ${VERSIONTAG}")
configure_file("src/version.cpp.in" "version.cpp")
add_library(version version.cpp)

# MOSEK example
add_executable(test_lo1 test/lo1.cc)
target_link_libraries(test_lo1 ${MOSEK_C_LIBRARY} ${MOSEK_CXX_LIBRARY})
fix_mosek_lib_path(test_lo1)

# Utility library
add_library(util src/util.cpp)
target_link_libraries(util geometry-central)

# Handle-cutting library
add_library(handlecut src/handlecut.cpp)
target_link_libraries(handlecut util)

# Surface co-cutting library
add_library(cocut src/cocut.cpp)
target_link_libraries(cocut util)

# Untangling library
add_library(untangle src/untangle.cpp)
target_link_libraries(untangle ultimaille)

# G-mapping library
add_library(gmap src/gmap.cpp)
target_link_libraries(gmap util ${MOSEK_C_LIBRARY} ${MOSEK_CXX_LIBRARY} untangle)

# Path lifting library
add_library(lift src/transformation.cpp src/lift.cpp)
target_link_libraries(lift util)

# Refinement library
add_library(refine src/refine.cpp)
target_link_libraries(refine util)

# Tutte embedding test
add_executable(test_tutte test/tutte.cpp)
target_link_libraries(test_tutte ${MOSEK_C_LIBRARY} ${MOSEK_CXX_LIBRARY})
fix_mosek_lib_path(test_tutte)

# Isometric mapping test
add_executable(test_isometric test/isometric.cpp)
target_link_libraries(test_isometric util ${MOSEK_C_LIBRARY} ${MOSEK_CXX_LIBRARY})
fix_mosek_lib_path(test_isometric)

# Conformal mapping test
add_executable(test_conformal test/conformal.cpp)
target_link_libraries(test_conformal util ${MOSEK_C_LIBRARY} ${MOSEK_CXX_LIBRARY})
fix_mosek_lib_path(test_conformal)

# G-flattening test
add_executable(test_gflatten test/gflatten.cpp)
target_link_libraries(test_gflatten util ${MOSEK_C_LIBRARY} ${MOSEK_CXX_LIBRARY})
fix_mosek_lib_path(test_gflatten)

# Surface cutting test
add_executable(test_cut test/cut.cpp)
target_link_libraries(test_cut util)

# Surface co-cutting test
add_executable(test_cocut test/cocut.cpp)
target_link_libraries(test_cocut cocut)

# G-mapping test
add_executable(test_gmap test/gmap.cpp)
target_link_libraries(test_gmap gmap)
fix_mosek_lib_path(test_gmap)

# Path lifting test
add_executable(test_lift test/lift.cpp)
target_link_libraries(test_lift lift)

# Refinement test
add_executable(test_refine test/refine.cpp)
target_link_libraries(test_refine refine)

# Handle-cutting test
add_executable(test_handlecut test/handlecut.cpp)
target_link_libraries(test_handlecut handlecut)

# Handle-cutter app
add_executable(app_handlecut src/app_handlecut.cpp)
target_link_libraries(app_handlecut handlecut version)
set_target_properties(app_handlecut PROPERTIES OUTPUT_NAME handlecut)   # https://gist.github.com/jlgerber/eafc4ee2b9954e27dd2bb009496b1b03

# Landmarker app
add_executable(landmarker src/app_landmarker.cpp)
target_link_libraries(landmarker glew imgui handlecut version)

# Seamless Surface Mapping app
add_executable(seamlesssurfmap src/app_seamlesssurfmap.cpp)
target_link_libraries(seamlesssurfmap cocut gmap lift refine version)
fix_mosek_lib_path(seamlesssurfmap)

# Untangle app
add_executable(app_untangle src/app_untangle.cpp)
target_link_libraries(app_untangle untangle util version)
set_target_properties(app_untangle PROPERTIES OUTPUT_NAME untangle)   # https://gist.github.com/jlgerber/eafc4ee2b9954e27dd2bb009496b1b03

# GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/1b18723e874b256c1e39378c6774a90701d70f7a.zip
)
FetchContent_MakeAvailable(googletest)
enable_testing()
include(GoogleTest)

# Unit test for util.cpp
add_executable(test_util ${CMAKE_CURRENT_LIST_DIR}/test/util.cpp)
target_link_libraries(test_util util gtest_main)
gtest_discover_tests(test_util)
