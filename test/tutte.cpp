#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>
#include <fusion.h>

#include "common.h"

void mosek_harmonic(const MatrixXd& V, const MatrixXi& F, const VectorXi& b, const VectorXd& bc, VectorXd& x_out) {
  MosekModel::t M = new MosekModel("harmonic"); auto _M = finally([&]() { M->dispose(); });

  const int nV = V.rows();
  const int nF = F.rows();

  Variable::t x = M->variable("x", nV, Domain::unbounded());
  Variable::t y = M->variable("y", nF * 3, Domain::unbounded());

  MatrixXd C;
  igl::cotmatrix_entries(V, F, C);

  // Energy
  for (int i = 0; i < nF; ++i) {
    for (int j = 0; j < 3; ++j) {
      int v0 = F(i, j);
      int v1 = F(i, (j + 1) % 3);
      double w = C(i, (j + 2) % 3);
      auto d = Expr::mul(w <= 0 ? 0 : std::sqrt(w), Expr::sub(x->index(v0), x->index(v1)));
      auto z = Expr::vstack(0.5, y->index(3 * i + j), d);
      M->constraint(std::string("qc") + std::to_string(3 * i + j), z, Domain::inRotatedQCone());
    }
  }
  M->objective("obj", ObjectiveSense::Minimize, Expr::sum(y));

  // Boundary condition
  for (int i = 0; i < b.size(); ++i) {
    M->constraint(std::string("bc") + std::to_string(i), x->index(b[i]), Domain::equalsTo(bc[i]));
  }

  // Solve & retrieve
  M->solve();

  // Get the linear solution values
  ndarray<double, 1> xlvl = *(x->level());

  x_out.resize(nV);
  for (int i = 0; i < nV; ++i)
    x_out[i] = xlvl[i];
}

int main(int argc, char *argv[])
{
  MatrixXd V, V_uv;
  MatrixXi F;
  // Load a mesh in OFF format
  igl::read_triangle_mesh(DATA_PATH "/camelhead.off", V, F);

  // Find the open boundary
  VectorXi bnd;
  igl::boundary_loop(F,bnd);

  // Map the boundary to a circle, preserving edge proportions
  MatrixXd bnd_uv;
  igl::map_vertices_to_circle(V,bnd,bnd_uv);

  // Harmonic parametrization for the internal vertices
  std::string postfix;
#if 1
  igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);
  postfix = "_libigl";
#else
  VectorXd V_u;
  VectorXd V_v;
  mosek_harmonic(V, F, bnd, bnd_uv.col(0), V_u);
  mosek_harmonic(V, F, bnd, bnd_uv.col(1), V_v);
  V_uv.resize(V.rows(), 2);
  V_uv << V_u, V_v;
  postfix = "_mosek";
#endif

  // Scale UV to make the texture more clear
  V_uv *= 5;

  // Output mesh
  V_uv.conservativeResize(V.rows(), 3);
  V_uv.col(2).setZero();
  igl::writeOBJ("tutte" + postfix + ".obj", V_uv, F);
}
