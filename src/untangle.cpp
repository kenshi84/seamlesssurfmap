// Adapted from https://github.com/ssloy/invertible-maps/blob/main/cpp/untangle2d.cpp
#include <ultimaille/all.h>

#include <igl/internal_angles.h>

#include "untangle.h"
#include "util.h"

using namespace UM;

namespace {

double triangle_area_2d(vec2 a, vec2 b, vec2 c) {
  return .5 * ((b.y-a.y)*(b.x+a.x) + (c.y-b.y)*(c.x+b.x) + (a.y-c.y)*(a.x+c.x));
}

double triangle_aspect_ratio_2d(vec2 a, vec2 b, vec2 c) {
  double l1 = (b-a).norm();
  double l2 = (c-b).norm();
  double l3 = (a-c).norm();
  double lmax = std::max(l1, std::max(l2, l3));
  return lmax*(l1+l2+l3)/(4.*std::sqrt(3.)*triangle_area_2d(a, b, c));
}

double chi(double eps, double det) {
  if (det>0)
    return (det + std::sqrt(eps*eps + det*det))*.5;
  return .5*eps*eps / (std::sqrt(eps*eps + det*det) - det);
}

double chi_deriv(double eps, double det) {
  return .5+det/(2.*std::sqrt(eps*eps + det*det));
}

}

struct Untangle2D {
  Untangle2D(Triangles &mesh) : m(mesh), X(m.nverts()*2), lock(m.points), ref_tri(m), J(m), K(m), det(m), area(m) {
    for (int t : facet_iter(m)) {
      area[t] = m.util.unsigned_area(t);
      vec2 A,B,C;
      m.util.project(t, A, B, C);

      double ar = triangle_aspect_ratio_2d(A, B, C);
      if (ar>10) { // if the aspect ratio is bad, assign an equilateral reference triangle
        double a = ((B-A).norm() + (C-B).norm() + (A-C).norm())/3.; // edge length is the average of the original triangle
        area[t] = sqrt(3.)/4.*a*a;
        A = {0., 0.};
        B = {a, 0.};
        C = {a/2., std::sqrt(3.)/2.*a};
      }

      mat<2,2> ST = {{B-A, C-A}};
      ref_tri[t] = mat<3,2>{{ {-1,-1},{1,0},{0,1} }}*ST.invert_transpose();
    }
  }

  void lock_boundary_verts() {
    SurfaceConnectivity fec(m);
    for (int v : vert_iter(m))
      lock[v] = fec.is_boundary_vert(v);
  }

  void evaluate_jacobian(const std::vector<double> &X) {
    detmin = std::numeric_limits<double>::max();
    ninverted = 0;
#pragma omp parallel for reduction(min:detmin) reduction(+:ninverted)
    for (int t=0; t<m.nfacets(); t++) {
      mat<2,2> &J = this->J[t];
      J = {};
      for (int i=0; i<3; i++)
        for (int d : range(2))
          J[d] += ref_tri[t][i]*X[2*m.vert(t,i) + d];
      this->K[t] = { {{ +J[1].y, -J[1].x }, { -J[0].y, +J[0].x }} };  // dual basis
      det[t] = J.det();
      detmin = std::min(detmin, det[t]);
      ninverted += (det[t]<=0);
    }
  }

  bool go() {
    std::vector<SpinLock> spin_locks(X.size());
    eps = 1;
    evaluate_jacobian(X);
    if (debug>0) std::cerr <<  "number of inverted elements: " << ninverted << std::endl;
    for (int iter=0; iter<maxiter; iter++) {
      if (debug>0) std::cerr << "iteration #" << iter << std::endl;
      const LBFGS_Optimizer::func_grad_eval func = [&](const std::vector<double>& X, double& F, std::vector<double>& G) {
        std::fill(G.begin(), G.end(), 0);
        F = 0;
        evaluate_jacobian(X);
//#pragma omp parallel for reduction(vec_double_plus:G) reduction(+:F)
#pragma omp parallel for reduction(+:F)
        for (int t=0; t<m.nfacets(); t++) {
          double c1 = chi(eps, det[t]);
          double c2 = chi_deriv(eps, det[t]);

          double f = (J[t][0]*J[t][0] + J[t][1]*J[t][1])/c1;
          double g = (1+det[t]*det[t])/c1;
          F += ((1-theta)*f + theta*g)*area[t];

          for (int dim : range(2)) {
            vec2 a = J[t][dim]; // tangent basis
            vec2 b = K[t][dim]; // dual basis
            vec2 dfda = (a*2. - b*f*c2)/c1;
            vec2 dgda = b*(2*det[t]-g*c2)/c1;
            for (int i=0; i<3; i++) {
              int v = m.vert(t,i);
              if (lock[v]) continue;
              spin_locks[v*2+dim].lock();
              G[v*2+dim] += ((dfda*(1.-theta) + dgda*theta)*ref_tri[t][i])*area[t];
              spin_locks[v*2+dim].unlock();
            }
          }
        }
//              double GG = std::sqrt(std::transform_reduce(G.begin(), G.end(), G.begin(), 0.));
//              double XX = std::sqrt(std::transform_reduce(X.begin(), X.end(), X.begin(), 0.));
//              std::cerr << F << " " << ninverted << " " << GG/std::max(1., XX) << std::endl;
      };

      double E_prev, E;
      std::vector<double> trash(X.size());
      func(X, E_prev, trash);

      LBFGS_Optimizer opt(func);
      opt.gtol = bfgs_threshold;
      opt.maxiter = bfgs_maxiter;
      opt.verbose = debug > 0;
      opt.run(X);

      func(X, E, trash);
      if (debug>0) std::cerr << "E: " << E << " eps: " << eps << " detmin: " << detmin << " ninv: " << ninverted << std::endl;

#if 0
      double sigma = std::max(1.-E/E_prev, 1e-1);
      if (detmin>=0)
        eps *= (1-sigma);
      else
        eps *= 1 - (sigma*std::sqrt(detmin*detmin + eps*eps))/(std::abs(detmin) + std::sqrt(detmin*detmin + eps*eps));

#else
      double sigma = std::max(1.-E/E_prev, 1e-1);
      double mu = (1-sigma)*chi(eps, detmin);
      if (detmin<mu)
        eps = std::max(1e-9, 2*std::sqrt(mu*(mu-detmin)));
      else eps = 1e-9;
#endif

      if  (detmin>0 && std::abs(E_prev - E)/E<converge_threshold) break;
    }
    return !ninverted;
  }

  ////////////////////////////////
  // Untangle2D state variables //
  ////////////////////////////////

  // optimization input parameters
  Triangles &m;           // the mesh to optimize
  double theta = 1./128.; // the energy is (1-theta)*(shape energy) + theta*(area energy)
  int maxiter = 10000;    // max number of outer iterations
  double bfgs_threshold = 1e-4;
  int bfgs_maxiter = 30000; // max number of inner iterations
  int debug = 1;          // verbose level
  double converge_threshold = 1e-5;

  // optimization state variables

  std::vector<double> X;     // current geometry
  PointAttribute<bool> lock; // currently lock = boundary vertices
  FacetAttribute<mat<3,2>> ref_tri;
  FacetAttribute<mat<2,2>> J; // per-tet Jacobian matrix = [[JX.x JX.y, JX.z], [JY.x, JY.y, JY.z], [JZ.x, JZ.y, JZ.z]]
  FacetAttribute<mat<2,2>> K; // per-tet dual basis: det J = dot J[i] * K[i]
  FacetAttribute<double> det; // per-tet determinant of the Jacobian matrix
  FacetAttribute<double> area; // reference area
  double eps;       // regularization parameter, depends on min(jacobian)

  double detmin;    // min(jacobian) over all tetrahedra
  int ninverted; // number of inverted tetrahedra
};

bool untangle(const MatrixXd& V, const MatrixXi& F, MatrixXd& V_uv, UntangleParam* param) {
  const int nV = V.rows();
  const int nF = F.rows();
  ASSERT(V_uv.rows() == nV);

  MatrixXd V_uv_orig = V_uv;    // To restore original when failed

  Triangles m;

  // Copy vertices
  m.points.resize(nV);
  for (int i = 0; i < V.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      m.points[i][j] = V(i, j);
    }
  }

  // Copy faces
  m.facets.resize(3 * nF);
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      m.facets[3 * i + j] = F(i, j);
    }
  }

  // These are used to undo the scaling we apply to the model
  AlignedBox2d bb;
  const double boxsize = 10.;

  // Scale the target domain for better numerical stability
  for (int v : vert_iter(m)) {
    bb.extend((Vector2d)V_uv.row(v).transpose());
  }
  double maxside = bb.diagonal().maxCoeff();
  for (int v : vert_iter(m)) {
    Vector2d uv = V_uv.row(v).transpose();
    uv = (uv - (bb.max() + bb.min()) / 2.) * boxsize / maxside + Vector2d(1, 1) * boxsize / 2.;
    V_uv.row(v) = uv.transpose();
  }

  // Scale the input geometry to have the same area as the target domain
  double target_area = 0;
  for (int t : facet_iter(m)) {
    vec2 a = kt84::vector_cast<2, vec2>((RowVector2d)V_uv.row(m.vert(t, 0)));
    vec2 b = kt84::vector_cast<2, vec2>((RowVector2d)V_uv.row(m.vert(t, 1)));
    vec2 c = kt84::vector_cast<2, vec2>((RowVector2d)V_uv.row(m.vert(t, 2)));
    target_area += triangle_area_2d(a, b, c);
  }
  um_assert(target_area > 0); // Ascertain mesh requirements

  double source_area = 0;
  for (int t : facet_iter(m)) {
    source_area += m.util.unsigned_area(t);
  }

  for (vec3 &p : m.points) {
    p *= std::sqrt(target_area / source_area);
  }

  // Configure optimizer
  Untangle2D opt(m);
  if (param) {
    opt.theta = param->theta;
    opt.maxiter = param->maxiter;
    opt.bfgs_threshold = param->bfgs_threshold;
    opt.bfgs_maxiter = param->bfgs_maxiter;
    opt.debug = param->debug;
    opt.converge_threshold = param->converge_threshold;
  }

  // Copy from V_uv to opt.X
  for (int v : vert_iter(m)) {
    for (int d : range(2)) {
      opt.X[2 * v + d] = V_uv(v, d);
    }
  }

  // Lock boundary vertices
  opt.lock_boundary_verts();

  SPDLOG_INFO("Start untangling...");
  auto t1 = std::chrono::high_resolution_clock::now();

  // Go!
  bool success = opt.go();

  auto t2 = std::chrono::high_resolution_clock::now();
  SPDLOG_INFO("  Done, elapsed time: {}s", std::chrono::duration<double>(t2 - t1).count());

  // If successful, copy from opt.X to V_uv
  if (success) {
    for (int v : vert_iter(m)) {
      for (int d : range(2)) {
        V_uv(v, d) = opt.X[2 * v + d];
      }
    }

    // Check angle sum for each interior vertex
    MatrixXd K;
    igl::internal_angles(V_uv, F, K);

    std::vector<double> vertexAngleSums(nV, 0.0);
    for (int i = 0; i < F.rows(); ++i) {
      for (int j = 0; j < 3; ++j) {
        vertexAngleSums[F(i, j)] += K(i, j);
      }
    }

    SurfaceConnectivity fec(m);
    for (int v : vert_iter(m)) {
      if (!fec.is_boundary_vert(v)) {
        double d = vertexAngleSums[v] - 2.0 * util::PI();
        if (d > util::PI()) {
          SPDLOG_ERROR("  Vertex {} angle sum is larger than 2*pi by {}!", v, d);
          success = false;
        }
      }
    }
  }

  if (!success) {
    SPDLOG_ERROR("  Untangling failed, restoring original UV");
    V_uv = V_uv_orig;
    return false;
  }

  // Restore scale
  for (int v : vert_iter(m)) {
    Vector2d uv = V_uv.row(v).transpose();
    uv = (uv - Vector2d(1, 1) * boxsize / 2.) / boxsize * maxside + (bb.min() + bb.max()) / 2.;
    V_uv.row(v) = uv.transpose();
  }

  return true;
}
