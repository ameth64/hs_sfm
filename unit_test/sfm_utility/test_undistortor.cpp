#include <gtest/gtest.h>
#include <iostream>

#include "hs_sfm/sfm_utility/undistortor.hpp"
#include "hs_sfm/sfm_utility/radial_distortor.hpp"
#include "hs_sfm/sfm_utility/decentering_distortor.hpp"

namespace
{

template <typename _Scalar>
class TestUndistortor
{
public:
  typedef _Scalar Scalar;
  typedef hs::sfm::Undistortor<Scalar> Undistortor;
  typedef typename Undistortor::IntrinsicParams IntrinsicParams;
  typedef typename Undistortor::ImagePoint ImagePoint;
  typedef hs::sfm::RadialDistortor<Scalar> RadialDistortor;
  typedef hs::sfm::DecenteringDistortor<Scalar> DecenteringDistortor;

  bool Test(const IntrinsicParams& intrinsic_params,
            const ImagePoint& point_undistort_true,
            Scalar precision) const
  {
    typedef typename IntrinsicParams::KMatrix KMatrix;

    ImagePoint point_distort;
    KMatrix K = intrinsic_params.GetKMatrix();
    KMatrix K_inv = K.inverse();
    point_distort = K_inv.block(0, 0, 2, 2) * point_undistort_true +
                    K_inv.block(0, 2, 2, 1);

    Scalar radial_dx, radial_dy, decenter_dx, decenter_dy;
    RadialDistortor radial_distortor;
    radial_distortor(intrinsic_params.k1(),
                     intrinsic_params.k2(),
                     intrinsic_params.k3(),
                     point_distort[0], point_distort[1],
                     radial_dx, radial_dy);
    DecenteringDistortor dicentering_distortor;
    dicentering_distortor(intrinsic_params.d1(), intrinsic_params.d2(),
                          point_distort[0], point_distort[1],
                          decenter_dx, decenter_dy);
    point_distort[0] += radial_dx + decenter_dx;
    point_distort[1] += radial_dy + decenter_dy;

    point_distort = K.block(0, 0, 2, 2) * point_distort + K.block(0, 2, 2, 1);

    Undistortor undistortor;
    ImagePoint point_undistort =
      undistortor.UndistortImagePoint(point_distort, intrinsic_params);

    return point_undistort.isApprox(point_undistort_true, precision);
  }
};

TEST(TestUndistortor, SimpleTest)
{
  typedef double Scalar;
  typedef TestUndistortor<Scalar> Tester;
  typedef Tester::ImagePoint ImagePoint;
  typedef Tester::IntrinsicParams IntrinsicParams;

  Tester tester;

  IntrinsicParams intrinsic_params(4880.97,
                                   0,
                                   3024.6,
                                   1993.09,
                                   1,
                                   -0.0993962,
                                   0.125977,
                                   -0.0454404);

  ImagePoint image_point;

  image_point << 5430.48,
                 3845.64;
  ASSERT_EQ(true, tester.Test(intrinsic_params, image_point, 1e-3));

  image_point << 3010.47,
                 1983.73;
  ASSERT_EQ(true, tester.Test(intrinsic_params, image_point, 1e-3));
}

}
