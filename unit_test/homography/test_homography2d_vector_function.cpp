#include <gtest/gtest.h>

#include "hs_sfm/homography/homography2d_vector_function.hpp"
#include "hs_sfm/homography/homography2d_synthetic_data_generator.hpp"

namespace
{

template <typename _Scalar>
class TestHomography2dVectorFunction
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::homography::Homography2DVectorFunction<Scalar>
          VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;

public:
  Err operator() (const VectorFunction& vector_function,
                  const XVector& x,
                  const YVector& true_y) const
  {
    YVector y;
    if (vector_function(x, y) != 0)
    {
      return -1;
    }

    if (y.isApprox(true_y))
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }
};

TEST(TestHomography2dVectorFunction, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef hs::sfm::homography::Homography2DSyntheticDataGenerator<
            Scalar, ImageDimension> DataGenerator;
  typedef hs::sfm::homography::Homography2DVectorFunction<Scalar>
          VectorFunction;
  typedef VectorFunction::XVector XVector;
  typedef VectorFunction::YVector YVector;
  typedef TestHomography2dVectorFunction<Scalar> Test;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 1;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 1;
  Scalar north_west_angle = 60;
  size_t number_of_keys = 5000;
  Scalar points_height = 25;

  DataGenerator data_generator(focal_length_in_metre,
                               number_of_strips,
                               ground_resolution,
                               image_width,
                               image_height,
                               pixel_size,
                               lateral_overlap_ratio,
                               longitudinal_overlap_ratio,
                               scene_max_height,
                               camera_height_stddev,
                               camera_planar_stddev,
                               camera_rotation_stddev,
                               north_west_angle,
                               number_of_keys,
                               points_height);

  VectorFunction vector_function;
  XVector x;
  YVector y;
  ASSERT_EQ(0, data_generator(vector_function, x, y));
  
  Test test;
  ASSERT_EQ(0, test(vector_function, x, y));
}

}