#include "transformation.h"

Matrix2d Transformation::getA() const {
  Matrix2d ret;
  ret <<
    a(), -b(),
    b(), a();
  return ret;
}

void Transformation::setA(const Matrix2d& A) {
  ASSERT(A(0, 0) == A(1, 1));
  ASSERT(A(0, 1) == -A(1, 0));
  a() = A(0, 0);
  b() = A(1, 0);
}

Vector2d Transformation::apply(const Vector2d& p) const {
  return getA() * p + t;
}

Transformation Transformation::inverse() const {
  Transformation r;
  /*
  for all p,
  r.A * (A * p + t) + r.t = p
  <-->
  (r.A * A) * p + r.A * t + r.t = p
  <-->
  r.A = inv(A)
  r.t = -r.A * t;
  */

  // Linear part
  r.setA(getA().inverse());

  // Translation part
  r.t = -r.getA() * t;

  return r;
}

Transformation Transformation::make(const Vector2d& p0, const Vector2d& p1, const Vector2d& q0, const Vector2d& q1) {
  /*
  A * p0 + t = q0,
  A * p1 + t = q1
  -->
  A * (p1 - p0) = q1 - q0
  -->
  dp := p1 - p0,
  dq := q1 - q0,
  [a -b] * [dp.x] = [dq.x]
  [b  a]   [dp.y]   [dq.y]
  -->
  [dp.x -dp.y] * [a] = [dq.x]
  [dp.y  dp.x]   [b]   [dq.y]
  */
  Vector2d dp = p1 - p0;
  Vector2d dq = q1 - q0;
  Matrix2d P;
  P <<
    dp.x(), -dp.y(),
    dp.y(),  dp.x();

  Transformation r;
  r.ab = P.inverse() * dq;
  r.t = q0 - r.getA() * p0;
  return r;
}
