#include <gtest/gtest.h>
#include <util.h>

TEST(util, polar_decomposition) {
  for (int i = 0; i < 100; ++i) {
    Matrix2d A = Matrix2d::Random();
    if (A.determinant() <= 0) continue;

    Matrix2d R, Y;
    util::polar_decomposition(A, R, Y);

    EXPECT_LT((R.transpose() * R - Matrix2d::Identity()).norm(), 1.0e-12);
    EXPECT_LT((A - R * Y).norm(), 1.0e-12);
  }
}

TEST(util, append_zero_col) {
  MatrixXd Min(4, 2);
  Min <<
    1,2,
    3,4,
    5,6,
    7,8;
  MatrixXd Mout = util::append_zero_col(Min);

  MatrixXd Mout2(4, 3);
  Mout2 <<
    1,2,0,
    3,4,0,
    5,6,0,
    7,8,0;

  EXPECT_EQ(Mout, Mout2);
}
