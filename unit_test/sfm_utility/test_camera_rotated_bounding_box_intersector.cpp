#include <iostream>
#include <cmath>

#include <gtest/gtest.h>

#include "hs_math/geometry/euler_angles.hpp"

#include "hs_sfm/sfm_utility/camera_rotated_bounding_box_intersector.hpp"

namespace
{

TEST(TestCameraRotatedBoundingBoxIntersector, SimpleTest)
{
  typedef double Scalar;
  typedef hs::sfm::CameraRotatedBoundingBoxIntersector<Scalar> Intersector;
  typedef Intersector::IntrinsicParams IntrinsicParams;
  typedef Intersector::ExtrinsicParams ExtrinsicParams;
  typedef Intersector::Point Point;
  typedef Intersector::RMatrix RMatrix;
  typedef hs::math::geometry::EulerAngles<Scalar> EulerAngles;

  IntrinsicParams intrinsic_params(7000.0, 0.0, 3000.0, 2000.0, 1.0);
  ExtrinsicParams extrinsic_params;
  size_t width = 6000;
  size_t height = 4000;
  Point center, min, max;
  RMatrix rmatrix;
  EulerAngles euler_angles;
  Intersector intersector;

  //包围盒顶点均在四棱锥外面，但包围盒与四棱锥仍相交的情况
  center << 50.5493,
            25.6043,
            375.0;
  min << -756.1627,
         -432.7434,
         350.0;
  max << 857.2612,
         483.9519,
         400.0;
  min -= center;
  max -= center;
  euler_angles[0] = -34.0 / 180.0 * M_PI;
  euler_angles[1] = -27.0 / 180.0 * M_PI;
  euler_angles[2] = -51.0 / 180.0 * M_PI;
  rmatrix = euler_angles.ToOrthoRotMat<1, 2, 3, 1>();
  rmatrix.transposeInPlace();
  ASSERT_EQ(true, intersector(intrinsic_params, extrinsic_params,
                              width, height, center, min, max, rmatrix));
  //包围盒顶点均在四棱锥里面的情况
  center << 3.5432,
            -4.49,
            365.0;
  min << -57.1568,
         -48.5698,
         350.0;
  max << 64.2433,
         39.5899,
         380.0;
  min -= center;
  max -= center;
  euler_angles[0] = -15.0 / 180.0 * M_PI;
  euler_angles[1] = -31.0 / 180.0 * M_PI;
  euler_angles[2] = -26.0 / 180.0 * M_PI;
  rmatrix = euler_angles.ToOrthoRotMat<1, 2, 3, 1>();
  rmatrix.transposeInPlace();
  ASSERT_EQ(true, intersector(intrinsic_params, extrinsic_params,
                              width, height, center, min, max, rmatrix));
  //包围盒部分顶点在四棱锥里面部分顶点在四棱锥外面的情况
  center << -52.7403,
            46.9207,
            372.5;
  min << -169.724,
         -48.5698,
         350.0;
  max << 64.2433,
         142.4111,
         395.0;
  min -= center;
  max -= center;
  euler_angles[0] = -15.0 / 180.0 * M_PI;
  euler_angles[1] = -31.0 / 180.0 * M_PI;
  euler_angles[2] = -26.0 / 180.0 * M_PI;
  rmatrix = euler_angles.ToOrthoRotMat<1, 2, 3, 1>();
  rmatrix.transposeInPlace();
  ASSERT_EQ(true, intersector(intrinsic_params, extrinsic_params,
                              width, height, center, min, max, rmatrix));
  //包围盒与四棱锥相离的情况
  center << 369.6993,
            -30.9792,
            365.0;
  min << 245.9582,
         -137.5126,
         350.0;
  max << 493.4403,
         75.5542,
         380.0;
  min -= center;
  max -= center;
  euler_angles[0] = -20.0 / 180.0 * M_PI;
  euler_angles[1] = -10.0 / 180.0 * M_PI;
  euler_angles[2] = -30.0 / 180.0 * M_PI;
  rmatrix = euler_angles.ToOrthoRotMat<1, 2, 3, 1>();
  rmatrix.transposeInPlace();
  ASSERT_EQ(false, intersector(intrinsic_params, extrinsic_params,
                              width, height, center, min, max, rmatrix));
#if 0
  Point corners[8] =
  {
      Point(min[0], min[1], min[2]),
      Point(min[0], min[1], max[2]),
      Point(min[0], max[1], min[2]),
      Point(min[0], max[1], max[2]),
      Point(max[0], min[1], min[2]),
      Point(max[0], min[1], max[2]),
      Point(max[0], max[1], min[2]),
      Point(max[0], max[1], max[2])
  };
  for (size_t i = 0; i < 8; i++)
  {
    Point rotated_corner = rmatrix.transpose() * corners[i] + center;
    std::cout<<"rotated_corner["<<i<<"]:\n"<<rotated_corner<<"\n";
  }
#endif
}

}
