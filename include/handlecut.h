#include "common.h"

enum handlecut_errorcode {
  HANDLECUT_SUCCESS = 0,
  HANDLECUT_PATH_NOTFOUND,
  HANDLECUT_LOOP_INVALID,
};

// Cuts along a loop going through the seed vertex.
// Returns true if successful, where the number of connected components remains the same, and the genus decreases by one.
bool handlecut(
  const MeshGeometry& input,
  int seedVertex,
  MeshGeometry& output,
  bool tunnelMode = false,
  int* correspondingVertex = nullptr,     // Newly created boundary vertex corresponding to seedVertex
  handlecut_errorcode* ec = nullptr);
