#include "common.h"
#include <fusion.h>

namespace mosek_util {

inline std::shared_ptr<ndarray<double, 2>> from_eigen_to_mosek(const Matrix2d& m) {
  return monty::new_array_ptr<double, 2>({
    {m(0, 0), m(0, 1)},
    {m(1, 0), m(1, 1)}
  });
}

inline std::shared_ptr<ndarray<double, 1>> from_eigen_to_mosek(const Vector2d& v) {
  return new_array_ptr<double>({ v[0], v[1] });
}

inline MosekMatrix::t get_similarity_matrix(const Vector2d& v) {
  double a = v[0];
  double b = v[1];
  return MosekMatrix::dense(monty::new_array_ptr<double, 2>({
    {a , -b},
    {b, a}
  }));
}

inline Expression::t get_variable_similarity_matrix(Variable::t v) {
  auto a = v->index(0);
  auto b = v->index(1);
  return Expr::vstack(
    Expr::hstack(a, Expr::neg(b)),
    Expr::hstack(b, a));
}

inline Expression::t trace(const Expression::t& m) {
  return Expr::sum(Expr::mulDiag(m, MosekMatrix::eye(2)));
}

}
