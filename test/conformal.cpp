#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>

#include <vector>

#include "util.h"
#include "mosek_util.h"

// Local per-triangle computation of map differential & best fit rotation
void conformal_local(const MatrixXd& V, const MatrixXi& F, const MatrixXd& V_uv, std::vector<util::MapDifferential>& mapDiff, std::vector<Matrix2d>& R, std::vector<double>& c) {
  const int nF = F.rows();

  mapDiff.resize(nF);
  R.resize(nF);
  c.resize(nF);

  for (int i = 0; i < nF; ++i) {
    mapDiff[i] = util::compute_map_differential(V, F, V_uv, i);

    Matrix2d Y;
    util::polar_decomposition(mapDiff[i].A, R[i], Y);

    double detA = mapDiff[i].A.determinant();
    ASSERT(detA > 0);
    c[i] = std::min<double>(100.0, std::max<double>(0.01, 1.0 / std::sqrt(detA)));
  }
}

// Global optimization via SOCP
void conformal_global(const MatrixXd& V, const MatrixXi& F, const std::vector<util::MapDifferential>& mapDiff, const std::vector<Matrix2d>& R, const std::vector<double>& c, MatrixXd& V_uv) {
  MosekModel::t M = new MosekModel("conformal");
  auto _M = finally([&]() { M->dispose(); });

  const int nV = V.rows();
  const int nF = F.rows();

  ASSERT(V.cols() == 3);
  ASSERT(F.cols() == 3);
  ASSERT(mapDiff.size() == nF);
  ASSERT(R.size() == nF);
  ASSERT(c.size() == nF);

  V_uv.resize(nV, 2);

  const double sqrt2 = std::sqrt(2.0);

  // Main variable (UV coordinate) per vertex
  Variable::t uv = M->variable("uv", 2 * nV, Domain::unbounded());

  // Auxiliary variable per triangle
  Variable::t rB = M->variable("rB", nF, Domain::unbounded());
  Variable::t rC = M->variable("rC", nF, Domain::unbounded());
  Variable::t u  = M->variable("u" , nF, Domain::unbounded());
  Variable::t s  = M->variable("s" , nF, Domain::unbounded());
  Variable::t S  = M->variable("S" , nF, Domain::unbounded());
  Variable::t U  = M->variable("U" , nF, Domain::unbounded());

  std::vector<double> sqrt_area(nF);

  SPDLOG_INFO("Building model...");

  for (int i = 0; i < nF; ++i) {
    auto rB_i = rB->index(i);
    auto rC_i = rC->index(i);
    auto u_i  = u ->index(i);
    auto s_i  = s ->index(i);
    auto S_i  = S ->index(i);
    auto U_i  = U ->index(i);

    auto Pinv_MOSEK = mosek_util::from_eigen_to_mosek(mapDiff[i].Pinv);

    sqrt_area[i] = std::sqrt(0.5 * mapDiff[i].P.determinant());

    // Variable 2D vectors representing the triangle corners
    int v0 = F(i, 0);
    int v1 = F(i, 1);
    int v2 = F(i, 2);

    auto q0 = Expr::vstack(uv->index(v0), uv->index(nV + v0));
    auto q1 = Expr::vstack(uv->index(v1), uv->index(nV + v1));
    auto q2 = Expr::vstack(uv->index(v2), uv->index(nV + v2));

    auto q01 = Expr::sub(q1, q0);
    auto q02 = Expr::sub(q2, q0);

    auto Q = Expr::hstack(q01, q02);

    auto A = Expr::mul(Q, Pinv_MOSEK);

    // Decompose A = B + C where B = 0.5 (A - A^T + tr(A)I)
    auto At = Expr::transpose(A);
    auto trA = mosek_util::trace(A);
    auto trA_I = Expr::mul(trA, MosekMatrix::eye(2));
    auto B = Expr::mul(0.5, Expr::sub(Expr::add(A, trA_I), At));
    auto C = Expr::sub(A, B);

    // Equation (10a) --------------+
    //    |B| + |C| <= sqrt(2) r    |
    //------------------------------+
    // |B| <= rB
    M->constraint("qc_10a_B_" + std::to_string(i),
      Expr::vstack(rB_i, Expr::flatten(B)),
      Domain::inQCone());

    // |C| <= rC
    M->constraint("qc_10a_C_" + std::to_string(i),
      Expr::vstack(rC_i, Expr::flatten(C)),
      Domain::inQCone());

    // r = (rB + rC) / sqrt(2)
    auto r = Expr::mul(1.0 / sqrt2, Expr::add(rB_i, rC_i));

    // Equation (10b) --------------------+
    //    |C| <= tr(R^T B) - sqrt(2) u    |
    //------------------------------------+
    auto R_MOSEK = MosekMatrix::dense(mosek_util::from_eigen_to_mosek(R[i]));
    auto RtB = Expr::mul(R_MOSEK->transpose(), B);
    auto tr_RtB = mosek_util::trace(RtB);

    M->constraint("qc_10b_" + std::to_string(i),
      Expr::vstack(Expr::sub(tr_RtB, Expr::mul(sqrt2, u_i)), Expr::flatten(C)),
      Domain::inQCone());

    // Equation (10c) --------------------+
    //    sqrt((u - s)^2 + 4) <= u + s    |
    //------------------------------------+
    M->constraint("qc_10c_" + std::to_string(i),
      Expr::vstack(Expr::add(u_i, s_i), Expr::sub(u_i, s_i), 2.0),
      Domain::inQCone());

    // Equation (10d) ------------------------+
    //    sqrt((s / c)^2 + (c * r)^2) <= S    |
    //----------------------------------------+
    M->constraint("qc_10d_" + std::to_string(i),
      Expr::vstack(S_i, Expr::mul(1.0 / c[i], s_i), Expr::mul(c[i], r)),
      Domain::inQCone());

    // Equation (10e) --------------------------+
    //    sqrt(S^2 + (U - 1/4)^2) <= U + 1/4    |
    //------------------------------------------+
    M->constraint("qc_10e_" + std::to_string(i),
      Expr::vstack(Expr::add(U_i, 0.25), S_i, Expr::sub(U_i, 0.25)),
      Domain::inQCone());
  }

  // Final SOCP constraint ---------------------+
  //    y >= a_1 * U_1^2 + ... + a_n * U_n^2    |
  //--------------------------------------------+
  Variable::t y = M->variable("y", Domain::unbounded());
  M->constraint("qc_y",
    Expr::vstack(0.5, y, Expr::mulElm(U, new_array_ptr(sqrt_area))),
    Domain::inRotatedQCone());

  // Fix one vertex to (0, 0)
  M->constraint("fix_uv_0",
    Expr::vstack(uv->index(0), uv->index(nV)),
    Domain::equalsTo(0.0));

  // Fix another vertex to lie on the line connecting vertex 0 & 1
  Vector2d dir = V_uv.row(1) - V_uv.row(0);
  auto q0 = Expr::vstack(uv->index(0), uv->index(nV));
  auto q1 = Expr::vstack(uv->index(1), uv->index(nV + 1));
  M->constraint("fix_uv_1",
    Expr::dot(Expr::sub(q1, q0), new_array_ptr<double>({ -dir[1], dir[0] })),
    Domain::equalsTo(0.0));

  // Set the objective function to minimize y
  M->objective("obj", ObjectiveSense::Minimize, Expr::mul(0.25, y));

  SPDLOG_INFO("Solving...");

  // Solve
  M->solve();
  SPDLOG_INFO("Primary objective value: " ANSI_BOLD ANSI_CYAN "{}" ANSI_RESET, M->primalObjValue());

  // Retrieve solution values
  ndarray<double, 1> uv_level = *(uv->level());
  for (int i = 0; i < nV; ++i)
    V_uv.row(i) << uv_level[i], uv_level[nV + i];
}

void tutte(const MatrixXd& V, const MatrixXi& F, MatrixXd& V_uv) {
  // Find the open boundary
  VectorXi bnd;
  igl::boundary_loop(F, bnd);

  // Map the boundary to a circle, preserving edge proportions
  MatrixXd bnd_uv;
  igl::map_vertices_to_circle(V, bnd, bnd_uv);

  // Harmonic parametrization for the internal vertices
  igl::harmonic(F, bnd, bnd_uv, 1, V_uv);
}

int main(int argc, char *argv[]) {
  // Load a mesh in OFF format
  MatrixXd V;
  MatrixXi F;
  igl::read_triangle_mesh(DATA_PATH "/bumpy-round-cube/A-disk.obj", V, F);

  // Initial state via Tutte embedding
  MatrixXd V_uv;
  tutte(V, F, V_uv);

  util::write_flattened_mesh(V_uv, F, "conformal-init.obj");

  SPDLOG_INFO("Initial state conformal energy: " ANSI_BOLD ANSI_YELLOW "{}" ANSI_RESET, util::compute_conformal_energy(V, F, V_uv));

  // Iterative optimization
  double prev_conformal_energy;
  kt84::MinSelector<int> bestIter;
  for (int iter = 0; ; ++iter) {
    SPDLOG_INFO("Iteration {}...", iter);

    std::vector<util::MapDifferential> mapDiff;
    std::vector<Matrix2d> R;
    std::vector<double> c;
    SPDLOG_INFO("Local step...");
    conformal_local(V, F, V_uv, mapDiff, R, c);

    if (iter > 0) {
      // Consistency check for conformal energy computation via isometric energy with optimal scaling
      double another_conformal_energy = util::compute_isometric_energy(V, F, V_uv, c);
      ASSERT(std::abs(prev_conformal_energy - another_conformal_energy) < 1.0e-10);
    }

    if (iter == 0) c = std::vector<double>(F.rows(), 1.0);     // Initial values should be 1

    SPDLOG_INFO("Global step...");
    conformal_global(V, F, mapDiff, R, c, V_uv);

    // Normalize total area
    double total_area = 0;
    for (int j = 0; j < F.rows(); ++j) {
      Vector2d q0 = V_uv.row(F(j, 0));
      Vector2d q1 = V_uv.row(F(j, 1));
      Vector2d q2 = V_uv.row(F(j, 2));
      Matrix2d Q;
      Q << q1 - q0, q2 - q0;
      total_area += 0.5 * Q.determinant();
    }
    V_uv *= 1.0 / std::sqrt(total_area);

    SPDLOG_INFO("Optimized conformal energy: " ANSI_BOLD ANSI_YELLOW "{}" ANSI_RESET, (prev_conformal_energy = util::compute_conformal_energy(V, F, V_uv)));

    util::write_flattened_mesh(V_uv, F, "conformal-" + std::to_string(iter) + ".obj");

    bestIter.update(prev_conformal_energy, iter);
    SPDLOG_INFO("Best iteration found so far: " ANSI_BOLD ANSI_MAGENTA "{}" ANSI_RESET, bestIter.value);
  }
}
