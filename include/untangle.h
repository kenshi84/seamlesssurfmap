#include "common.h"

struct UntangleParam {
  double theta = 1./128.; // the energy is (1-theta)*(shape energy) + theta*(area energy)
  int maxiter = 10000;    // max number of outer iterations
  double bfgs_threshold = 1e-4;
  int bfgs_maxiter = 30000; // max number of inner iterations
  int debug = 1;          // verbose level
  double converge_threshold = 1e-5;
};

bool untangle(const MatrixXd& V, const MatrixXi& F, MatrixXd& V_uv, UntangleParam* param = nullptr);
