#include <gtest/gtest.h>

#include "hs_sfm/bundle_adjustment/ba_naive_synthetic_data_generator.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"

namespace
{

TEST(TestBANaiveVectorFunction, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef hs::sfm::ba::BANaiveSyntheticDataGenerator<Scalar, ImageDimension>
          DataGenerator;

  typedef hs::sfm::ba::BANaiveVectorFunction<Scalar> VectorFunction;
  typedef VectorFunction::Index Index;
  typedef VectorFunction::XVector XVector;
  typedef VectorFunction::YVector YVector;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 10;
  size_t number_of_cameras_in_strip = 10;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 2000;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 1;
  Scalar north_west_angle = 60;

  VectorFunction vector_function;
  XVector x;
  YVector y, y1;

  DataGenerator data_generator(focal_length_in_metre,
                               number_of_strips,
                               number_of_cameras_in_strip,
                               ground_resolution,
                               image_width,
                               image_height,
                               pixel_size,
                               number_of_points,
                               lateral_overlap_ratio,
                               longitudinal_overlap_ratio,
                               scene_max_height,
                               camera_height_stddev,
                               camera_planar_stddev,
                               camera_rotation_stddev,
                               north_west_angle);
  ASSERT_EQ(0, data_generator(vector_function, x, y));

  ASSERT_EQ(0, vector_function(x, y1));

  for (Index i = 0; i < y1.rows(); i++)
  {
    ASSERT_NEAR(y[i], y1[i], Scalar(1e-5));
  }
}

}
