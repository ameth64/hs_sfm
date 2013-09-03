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
  typedef typename VectorFunction::XVec XVector;
  typedef typename VectorFunction::YVec YVector;
  
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
  typedef TestMultipleViewVectorFunction<Scalar> Test;
  typedef Test::VectorFunction VectorFunction;
  typedef VectorFunction::XVec XVec;
  typedef VectorFunction::YVec YVec;

  Scalar f = 0.019;
  size_t strip_num = 10;
  size_t cams_num_in_strip = 10;
  Scalar ground_resolution = 0.1;
  ImgDim img_width = 6000;
  ImgDim img_height = 4000;
  Scalar pixel_size = 0.0000039;
  Scalar lateral_overlap = 0.6;
  Scalar longitudinal_overlap = 0.8;
  Scalar scene_max_height = 50;
  Scalar cam_height_stddev = 5;
  Scalar cam_plannar_stddev = 5;
  Scalar cam_rot_stddev = 1;
  Scalar nw_angle = 60;

  VectorFunction vector_function;
  XVec x;
  YVec y;

  DataGenerator data_generator(f, strip_num, cams_num_in_strip,
                               ground_resolution, img_width, img_height,
                               pixel_size, lateral_overlap,
                               longitudinal_overlap,
                               scene_max_height,
                               cam_height_stddev,
                               cam_plannar_stddev,
                               cam_rot_stddev,
                               nw_angle);
  ASSERT_EQ(0, data_generator(vector_function, x, y));
  ASSERT_EQ(0, Test::Test(vector_function, x, y));
          
}
}