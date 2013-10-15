#include <iostream>

#include <gtest/gtest.h>

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_synthetic_data_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_relative_data_generator.hpp"

namespace
{

template <typename _Scalar>
class TestBAGCPConstrainedRelativeDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::ba::BAGCPConstrainedRelativeDataGenerator<Scalar> Generator;
  typedef hs::sfm::ba::BAGCPConstrainedVectorFunction<Scalar> GCPVectorFunction;
  typedef hs::sfm::ba::BANaiveVectorFunction<Scalar> NaiveVectorFunction;
  typedef typename Generator::XVector XVector;
  typedef typename Generator::YVector YVector;
  typedef typename Generator::Index Index;
  typedef typename Generator::Rotation Rotation;
  typedef typename Generator::Translate Translate;

public:
  Err operator()(const GCPVectorFunction& gcp_vector_function,
                 const XVector& absolute_x) const
  {
    NaiveVectorFunction naive_vector_function(gcp_vector_function);
    YVector absolute_naive_y;
    if (naive_vector_function(absolute_x, absolute_naive_y) != 0)
    {
      std::cout<<"naive vector function failed!\n";
      return -1;
    }

    Scalar scale;
    Rotation rotation;
    Translate translate;

    Generator generator;
    XVector relative_x;
    if (generator(gcp_vector_function, absolute_x,
                  relative_x,
                  scale,
                  rotation,
                  translate) != 0)
    {
      std::cout<<"generator failed!\n";
    }

    YVector relative_naive_y;
    if (naive_vector_function(relative_x, relative_naive_y) != 0)
    {
      std::cout<<"naive vector function failed!\n";
      return -1;
    }

    const Scalar threshold = Scalar(1e-8);
    Err result = 0;
    YVector diff_y = relative_naive_y - absolute_naive_y;
    for (Index i = 0; i < relative_naive_y.rows(); i++) 
    {
      Scalar relative_value = relative_naive_y[i];
      Scalar absolute_value = absolute_naive_y[i];
      Scalar diff = std::abs(diff_y[i]);
      if (diff > threshold)
      {
        std::cout<<"difference too large!";
        std::cout<<"relative_naive_y["<<i<<"] = "<<relative_value<<"\n";
        std::cout<<"absolute_naive_y["<<i<<"] = "<<absolute_value<<"\n";
        result = -1;
      }
    }

    return result;
  }
};

TEST(TestBAGCPConstrainedRelativeDataGenerator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef
    hs::sfm::ba::BAGCPConstrainedSyntheticDataGenerator<Scalar, ImageDimension>
    DataGenerator;
  typedef TestBAGCPConstrainedRelativeDataGenerator<Scalar> Test;
  typedef EIGEN_MATRIX(Scalar, 2, 2) Matrix22;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef Test::GCPVectorFunction GCPVectorFunction;
  typedef GCPVectorFunction::XVector XVector;
  typedef GCPVectorFunction::YVector YVector;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 3;
  size_t number_of_cameras_in_strip = 3;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 30;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 1;
  Scalar north_west_angle = 60;
  size_t number_of_gcps = 5;

  GCPVectorFunction gcp_vector_function;
  XVector x;
  YVector y;

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
                               north_west_angle,
                               number_of_gcps);

  ASSERT_EQ(0, data_generator(gcp_vector_function, x, y));

  Test test;
  ASSERT_EQ(0, test(gcp_vector_function, x));
}

}