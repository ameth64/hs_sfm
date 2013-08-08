#include <gtest/gtest.h>

#include "hs_sfm/utility/cam_type.hpp"

namespace
{

TEST(TestSFMType, SimpleTest)
{
  typedef double Scalar;

  typedef hs::sfm::IntrinParam<Scalar> Intrinc;

  Intrinc i0;
  Intrinc i1(2000);
  Intrinc i2(2000, 0.5, 20, 21, 1.5);
  Intrinc::KMat K0 = i2.getKMat();
  Intrinc i3(K0);
  Intrinc::KMat K1 = Intrinc::KMat(i3);
  Intrinc i4(i3);
  Intrinc i5 = i4;
}

}
