#include <iostream>
#include <gtest/gtest.h>

#include "hs_sfm/homography/homography2d_dlt_calculator.hpp"

TEST(TestHomography2DDLTCalculator, ExactTest)
{
  typedef double Scalar;
  typedef hs::sfm::homography::Homography2DDLTCalculator<Scalar> Calculator;
  typedef Calculator::Key Key;
  typedef Calculator::Correspondence Correspondence;
  typedef Calculator::CorrespondenceContainer CorrespondenceContainer;
  typedef Calculator::HMatrix HMatrix;
  typedef EIGEN_VECTOR(Scalar, 3) HKey;

  CorrespondenceContainer correspondences(4);
  correspondences[0].first << Scalar(-18),
                              Scalar(-13.5);
  correspondences[0].second << Scalar(-2800),
                               Scalar(-1950);
  correspondences[1].first << Scalar(18),
                              Scalar(-13.5);
  correspondences[1].second << Scalar(2750),
                               Scalar(-1847);
  correspondences[2].first << Scalar(18),
                              Scalar(13.5);
  correspondences[2].second << Scalar(2894),
                               Scalar(1977);
  correspondences[3].first << Scalar(-18),
                              Scalar(13.5);
  correspondences[3].second << Scalar(-2933),
                               Scalar(1954);

  HMatrix h_matrix;
  Calculator calculator;
  ASSERT_EQ(0, calculator(correspondences, h_matrix));
  for (size_t i = 0; i < correspondences.size(); i++)
  {
    HKey hkey_1;
    hkey_1.template segment<2>(0) = correspondences[i].first;
    hkey_1[2] = 1;
    HKey hkey_2;
    hkey_2.template segment<2>(0) = correspondences[i].second;
    hkey_2[2] = 1;
    HKey hkey_transformed = h_matrix * hkey_1;
    hkey_transformed /= hkey_transformed[2];
    ASSERT_EQ(true, hkey_transformed.isApprox(hkey_2, Scalar(1e-6)));
  }
}
