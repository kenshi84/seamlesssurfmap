#include "common.h"

struct Transformation {
  Vector2d ab = {0, 0};     // Linear part A = [a, -b; b, a]
  Vector2d t = {0, 0};      // Translation part [tx; ty]

  double& a() { return ab[0]; }
  double& b() { return ab[1]; }
  const double& a() const { return ab[0]; }
  const double& b() const { return ab[1]; }

  Matrix2d getA() const;

  void setA(const Matrix2d& A);

  Vector2d apply(const Vector2d& p) const;

  Transformation inverse() const;

  static Transformation make(const Vector2d& p0, const Vector2d& p1, const Vector2d& q0, const Vector2d& q1);
};

