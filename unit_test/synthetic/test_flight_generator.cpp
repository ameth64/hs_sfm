#include <iostream>

#include "gtest/gtest.h"

#include "hs_sfm/synthetic/flight_generator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestFlightGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;
  typedef hs::sfm::synthetic::FlightGenerator<Scalar, ImageDimension>
          Generator;

private:
  typedef typename Generator::ExtrinsicParams ExtrinsicParams;
  typedef typename Generator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename Generator::Point3D Point3D;
  typedef typename Generator::Point3DContainer Point3DContainer;
  typedef typename Generator::Image Image;
  typedef typename Generator::ImageContainer ImageContainer;

public:
  TestFlightGenerator(
    Scalar focal_length_in_metre,
    size_t number_of_strips,
    size_t number_of_cameras_in_strip,
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    size_t number_of_points,
    Scalar lateral_overlap_ratio,
    Scalar longitudinal_overlap_ratio,
    Scalar scene_max_height,
    Scalar camera_height_stddev,
    Scalar camera_planar_stddev,
    Scalar camera_rotation_stddev)
    : generator_(focal_length_in_metre,
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
                 camera_rotation_stddev) {}

public:
  Err Test() const
  {
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    Point3DContainer points;

    if (TestGenerator(extrinsic_params_set, images, points) != 0)
    {
      std::cout<<"test generator failed!\n";
      return -1;
    }

    if (TestNumberOfCameras(extrinsic_params_set,
                            images) != 0)
    {
      std::cout<<"test number of cameras failed!\n";
      return -1;
    }

    return 0;
  }

private:
  Err TestGenerator(ExtrinsicParamsContainer& extrinsic_params_set,
                    ImageContainer& images,
                    Point3DContainer& points) const
  {
    if (generator_(extrinsic_params_set,
                   images,
                   points) != 0)
    {
      return -1;
    }

    return 0;
  }

  Err TestNumberOfCameras(const ExtrinsicParamsContainer& extrinsic_params_set,
                          const ImageContainer& images) const
  {
    if (extrinsic_params_set.size() != images.size() ||
        extrinsic_params_set.size() !=
        generator_.number_of_strips() *
        generator_.number_of_cameras_in_strip())
    {
      return -1;
    }

    return 0;
  }

private:
  Generator generator_;
};

TEST(TestFlightGenerator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestFlightGenerator<Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 15;
  size_t number_of_cameras_in_strip = 20;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 2000;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 2;
  Scalar camera_planar_stddev = 2;
  Scalar camera_rotation_stddev = 10;

  Test test(focal_length_in_metre,
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
            camera_rotation_stddev);

  ASSERT_EQ(0, test.Test());
}

}
