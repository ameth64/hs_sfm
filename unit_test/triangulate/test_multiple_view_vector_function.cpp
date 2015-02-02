#include <iostream>

#include <gtest/gtest.h>

#include "test_multiple_view_base.hpp"

namespace
{
template <typename _Scalar>
class TestMultipleViewVectorFunction
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::triangulate::MultipleViewVectorFunction<Scalar>
          VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;

  static Err Test(const VectorFunction& vector_function,
                  const XVector& x,
                  const YVector& true_y)
  {
    YVector y;
    if (vector_function(x, y) != 0)
    {
      std::cout<<"vector function failed!\n";
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

TEST(TestMultipleViewVectorFunction, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImgDim;
  typedef hs::sfm::triangulate::MultipleViewSyntheicDataGenerator<Scalar,
                                                                  ImgDim>
          DataGenerator;
  typedef DataGenerator::IntrinsicParams IntrinsicParams;
  typedef DataGenerator::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef TestMultipleViewVectorFunction<Scalar> Test;
  typedef Test::VectorFunction VectorFunction;
  typedef VectorFunction::XVector XVector;
  typedef VectorFunction::YVector YVector;

  Scalar focal_length_in_metre = 0.035;
  size_t number_of_strips = 2;
  size_t number_of_cameras_in_strip = 3;
  Scalar ground_resolution = 0.1;
  ImgDim image_width = 6000;
  ImgDim image_height = 4000;
  Scalar pixel_size = 0.00000444313;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 1;
  Scalar north_west_angle = 60;

  IntrinsicParamsContainer intrinsic_params_set;
  IntrinsicParams intrinsic_params0(8506.46,
                                    0,
                                    2704.86,
                                    1762.62,
                                    1,
                                    0.0599996,
                                    -0.23131,
                                    -0.00460478,
                                    0.000162892,
                                    -0.000169231);
  intrinsic_params_set.push_back(intrinsic_params0);
  intrinsic_params_set.push_back(intrinsic_params0);
  intrinsic_params_set.push_back(intrinsic_params0);
  IntrinsicParams intrinsic_params1(8501.47,
                                    0,
                                    2752.73,
                                    1862.18,
                                    1,
                                    0.0576605,
                                    -0.254137,
                                    0.14565,
                                    -0.00114618,
                                    -0.000143283);
  intrinsic_params_set.push_back(intrinsic_params1);
  intrinsic_params_set.push_back(intrinsic_params1);
  intrinsic_params_set.push_back(intrinsic_params1);

  VectorFunction vector_function;
  XVector x;
  YVector y;

  DataGenerator data_generator(focal_length_in_metre,
                               number_of_strips,
                               number_of_cameras_in_strip,
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
                               intrinsic_params_set);
  ASSERT_EQ(0, data_generator(vector_function, x, y));
  ASSERT_EQ(0, Test::Test(vector_function, x, y));
}

}
