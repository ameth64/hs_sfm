#include <fstream>

#include <gtest/gtest.h>

#include <cereal/archives/portable_binary.hpp>

#include "hs_test_utility/test_toolkit/test_serialization.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"

namespace
{

TEST(TestCameraType, SimpleTest)
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

TEST(TestCameraType, SerializationTest)
{
  typedef double Scalar;
  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;

  IntrinsicParams intrinsic_params(4880.22,
                                   0.0,
                                   3026.11,
                                   1992.36,
                                   1.0,
                                   -0.100423,
                                   0.128487,
                                   -0.0482081,
                                   -3.17902e-05,
                                   7.75629e-05);

  hs::test::TestLocalSerialization<IntrinsicParams> intrinsic_params_tester;
  ASSERT_EQ(true, intrinsic_params_tester(intrinsic_params,
                                          "intrinsic_params.bin"));

  ExtrinsicParams extrinsic_params;
  extrinsic_params.rotation()[0] = 0.3;
  extrinsic_params.rotation()[1] = 0.4;
  extrinsic_params.rotation()[2] = 0.5;
  extrinsic_params.position()[0] = 0.7;
  extrinsic_params.position()[1] = 0.8;
  extrinsic_params.position()[2] = 0.9;

  hs::test::TestLocalSerialization<ExtrinsicParams> extrinsic_params_tester;
  ASSERT_EQ(true, extrinsic_params_tester(extrinsic_params,
                                          "extrinsic_params.bin"));

}

}
