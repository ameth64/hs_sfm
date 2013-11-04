#include <iostream>

#include <gtest/gtest.h>

#include "hs_sfm/projective/mle/intrinsic_constrained_synthetic_generator.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_vector_function.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestIntrinsicConstrainedVectorFunction
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
  typedef hs::sfm::projective::IntrinsicConstrainedSyntheticGenerator<
            Scalar, ImageDimension> SyntheticGenerator;
  typedef hs::sfm::projective::IntrinsicConstrainedVectorFunction<Scalar>
          VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;

public:
  TestIntrinsicConstrainedVectorFunction(
    Scalar focal_length_in_metre,
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    size_t number_of_points,
    Scalar scene_max_height,
    Scalar camera_height_stddev,
    Scalar camera_planar_stddev,
    Scalar camera_rot_stddev,
    Scalar north_west_angle,
    Scalar skew,
    Scalar principal_point_x,
    Scalar principal_point_y,
    Scalar pixel_ratio)
    : synthetic_generator_(focal_length_in_metre,
                           ground_resolution,
                           image_width,
                           image_height,
                           pixel_size,
                           number_of_points,
                           scene_max_height,
                           camera_height_stddev,
                           camera_planar_stddev,
                           camera_rot_stddev,
                           north_west_angle,
                           skew,
                           principal_point_x,
                           principal_point_y,
                           pixel_ratio) {}

  Err Test()
  {
    VectorFunction vector_function;
    XVector x;
    YVector y;
    if (synthetic_generator_(vector_function, x, y) != 0)
    {
      return -1;
    }

    YVector y_test;
    if (vector_function(x, y_test) != 0)
    {
      return -1;
    }

    if (!y.isApprox(y_test, 1e-10))
    {
      return -1;
    }

    return 0;
  }

private:
  SyntheticGenerator synthetic_generator_;
};

TEST(TestIntrinsicConstrainedVectorFunction, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestIntrinsicConstrainedVectorFunction<Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 5000;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rot_stddev = 1;
  Scalar north_west_angle = 60;
  Scalar skew = 0.0001;
  Scalar principal_point_x = 20.0;
  Scalar principal_point_y = 30.0;
  Scalar pixel_ratio = 1.0001;

  Test test(focal_length_in_metre,
            ground_resolution,
            image_width,
            image_height,
            pixel_size,
            number_of_points,
            scene_max_height,
            camera_height_stddev,
            camera_planar_stddev,
            camera_rot_stddev,
            north_west_angle,
            skew,
            principal_point_x,
            principal_point_y,
            pixel_ratio);

  ASSERT_EQ(0, test.Test());
}

}
