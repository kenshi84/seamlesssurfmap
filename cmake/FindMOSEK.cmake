# Adapted from https://github.com/libigl/libigl/blob/v2.3.0/cmake/FindMOSEK.cmake
#
# Try to find MOSEK
# Once done this will define
#
# MOSEK_FOUND           - system has MOSEK
# MOSEK_INCLUDE_DIR     - the MOSEK include directory
# MOSEK_C_LIBRARY       - MOSEK C API library
# MOSEK_CXX_LIBRARY     - MOSEK C++ Fusion API library
#

# Hardcoded search paths
set(SEARCH_PATHS
  ${CMAKE_SOURCE_DIR}/mosek/9.2/tools/platform/osx64x86/ 
  /usr/local/mosek/7/tools/platform/osx64x86/
  /usr/local/mosek/8/tools/platform/osx64x86/
  /usr/local/mosek/9.2/tools/platform/osx64x86/
  /usr/local/mosek/9.3/tools/platform/osx64x86/
  /usr/local/mosek/10.0/tools/platform/osx64x86/
  /usr/local/mosek/10.0/tools/platform/osxaarch64/
  /opt/mosek/7/tools/platform/linux64x86/
)

find_path(MOSEK_INCLUDE_DIR mosek.h
  PATHS ${SEARCH_PATHS}
  PATH_SUFFIXES h
)

find_library(MOSEK_C_LIBRARY NAMES mosek64
  HINT
    "${MOSEK_INCLUDE_DIR}"
    "${MOSEK_INCLUDE_DIR}/../bin"
    "${MOSEK_INCLUDE_DIR}/lib"
  PATHS
    ${SEARCH_PATHS}
  NO_DEFAULT_PATH
  PATH_SUFFIXES a bin lib dylib)

find_library(MOSEK_CXX_LIBRARY NAMES fusion64
  HINT
    "${MOSEK_INCLUDE_DIR}"
    "${MOSEK_INCLUDE_DIR}/../bin"
    "${MOSEK_INCLUDE_DIR}/lib"
  PATHS
    ${SEARCH_PATHS}
  NO_DEFAULT_PATH
  PATH_SUFFIXES a bin lib dylib)

# Check that Mosek was successfully found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  MOSEK DEFAULT_MSG MOSEK_C_LIBRARY MOSEK_CXX_LIBRARY MOSEK_INCLUDE_DIR
)

# Hide variables from CMake-Gui options
mark_as_advanced(MOSEK_C_LIBRARY MOSEK_CXX_LIBRARY MOSEK_INCLUDE_DIR)
