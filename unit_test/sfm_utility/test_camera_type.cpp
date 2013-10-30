#include <gtest/gtest.h>

#include "hs_sfm/sfm_utility/camera_type.hpp"

namespace
{

TEST(TestSFMType, SimpleTest)
{
  typedef double Scalar;

  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;

  IntrinsicParams i0;
  IntrinsicParams i1(2000);
  IntrinsicParams i2(2000, 0.5, 20, 21, 1.5);
  IntrinsicParams::KMatrix K0 = i2.GetKMatrix();
  IntrinsicParams i3(K0);
  IntrinsicParams::KMatrix K1 = IntrinsicParams::KMatrix(i3);
  IntrinsicParams i4(i3);
  IntrinsicParams i5 = i4;
}

}
