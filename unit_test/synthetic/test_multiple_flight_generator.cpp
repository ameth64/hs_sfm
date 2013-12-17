#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "gtest/gtest.h"

#include "hs_sfm/synthetic/multiple_flight_generator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestMultipleFlightGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;
  typedef hs::sfm::synthetic::MultipleFlightGenerator<Scalar, ImageDimension>
          Generator;

public:
  typedef typename Generator::ExtrinsicParams ExtrinsicParams;
  typedef typename Generator::ExtrinsicParamsSet ExtrinsicParamsSet;
  typedef typename Generator::ExtrinsicParamsSetContainer
                   ExtrinsicParamsSetContainer;
  typedef typename Generator::Image Image;
  typedef typename Generator::ImageSet ImageSet;
  typedef typename Generator::ImageSetContainer ImageSetContainer;
  typedef typename Generator::Point3D Point3D;
  typedef typename Generator::Point3DContainer Point3DContainer;
  typedef typename Generator::FlightGenerator FlightGenerator;
  typedef typename Generator::FlightGeneratorContainer
                   FlightGeneratorContainer;

public:
  TestMultipleFlightGenerator(
    Scalar flight_longitudinal_overlap_ratio,
    Scalar flight_lateral_overlap_ratio,
    Scalar north_west_angle,
    Scalar north_west_angle_stddev,
    Scalar offset_stddev,
    const FlightGeneratorContainer& flight_generators,
    size_t number_of_points,
    const std::string& test_name)
    : generator_(flight_longitudinal_overlap_ratio,
                 flight_lateral_overlap_ratio,
                 north_west_angle,
                 north_west_angle_stddev,
                 offset_stddev,
                 flight_generators),
      number_of_points_(number_of_points),
      test_name_(test_name) {}

  Err Test() const
  {
    Point3DContainer points;
    ExtrinsicParamsSetContainer extrinsic_params_sets;
    ImageSetContainer image_sets;
    if (TestGenerator(points, extrinsic_params_sets, image_sets) != 0)
    {
      std::cout<<"Testing generator failed!\n";
      return -1;
    }

    if (TestExtrinsicParamsSets(extrinsic_params_sets, image_sets) != 0)
    {
      std::cout<<"Testing extrinsic params sets failed!\n";
      return -1;
    }

    if (TestPoints(points) != 0)
    {
      std::cout<<"Testing points failed!\n";
      return -1;
    }

    return 0;
  }

private:
  Err TestGenerator(Point3DContainer& points,
                    ExtrinsicParamsSetContainer& extrinsic_params_sets,
                    ImageSetContainer& image_sets) const
  {
    if (generator_.GeneratePoints(number_of_points_, points) != 0) return -1;
    if (generator_.GenerateExtrinsicParamsSets(extrinsic_params_sets,
                                               image_sets) != 0) return -1;
    return 0;
  }

  Err TestExtrinsicParamsSets(
    const ExtrinsicParamsSetContainer& extrinsic_params_sets,
    const ImageSetContainer& image_sets) const
  {
    size_t number_of_image_sets = image_sets.size();
    if (number_of_image_sets != extrinsic_params_sets.size())
    {
      return -1;
    }
    for (size_t i = 0; i < number_of_image_sets; i++)
    {
      size_t number_of_images = image_sets[i].size();
      if (number_of_images != extrinsic_params_sets[i].size())
      {
        return -1;
      }

      std::stringstream ss;
      ss<<test_name_<<"_extrinsic_params_"<<i<<".txt";
      std::string path;
      ss>>path;

      std::ofstream file(path.c_str(), std::ios::out);
      if (!file.is_open())
      {
        return -1;
      }

      for (size_t j = 0; j < number_of_images; j++)
      {
        file<<j<<" "
            <<extrinsic_params_sets[i][j].position()[0]<<" "
            <<extrinsic_params_sets[i][j].position()[1]<<" "
            <<extrinsic_params_sets[i][j].position()[2]<<"\n";
      }
    }

    return 0;
  }

  Err TestPoints(const Point3DContainer& points) const
  {
    if (points.size() != number_of_points_)
    {
      return -1;
    }

    std::string path = test_name_ + "_points.txt";
    std::ofstream file(path.c_str(), std::ios::out);
    if (!file.is_open())
    {
      return -1;
    }
    for (size_t i = 0; i < number_of_points_; i++)
    {
      file<<i<<" "
          <<points[i][0]<<" "
          <<points[i][1]<<" "
          <<points[i][2]<<"\n";
    }

    return 0;
  }

private:
  Generator generator_;
  size_t number_of_points_;
  std::string test_name_;
};

TEST(TestMultipleFlightGenerator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestMultipleFlightGenerator<Scalar, ImageDimension> Test;
  typedef Test::FlightGenerator FlightGenerator;
  typedef Test::FlightGeneratorContainer FlightGeneratorContainer;

  FlightGeneratorContainer flight_generators;

  Scalar focal_length_in_metre_0 = 0.019;
  size_t number_of_strips_0 = 15;
  size_t number_of_cameras_in_strip_0 = 20;
  Scalar ground_resolution_0 = 0.1;
  ImageDimension image_width_0 = 6000;
  ImageDimension image_height_0 = 4000;
  Scalar pixel_size_0 = 0.0000039;
  size_t number_of_points_0 = 2000;
  Scalar lateral_overlap_ratio_0 = 0.6;
  Scalar longitudinal_overlap_ratio_0 = 0.8;
  Scalar scene_max_height_0 = 50;
  Scalar camera_height_stddev_0 = 2;
  Scalar camera_planar_stddev_0 = 2;
  Scalar camera_rotation_stddev_0 = 10;

  FlightGenerator flight_generator_0(
    focal_length_in_metre_0,
    number_of_strips_0,
    number_of_cameras_in_strip_0,
    ground_resolution_0,
    image_width_0,
    image_height_0,
    pixel_size_0,
    number_of_points_0,
    lateral_overlap_ratio_0,
    longitudinal_overlap_ratio_0,
    scene_max_height_0,
    camera_height_stddev_0,
    camera_planar_stddev_0,
    camera_rotation_stddev_0);
  flight_generators.push_back(flight_generator_0);

  Scalar focal_length_in_metre_1 = 0.030;
  size_t number_of_strips_1 = 20;
  size_t number_of_cameras_in_strip_1 = 15;
  Scalar ground_resolution_1 = 0.1;
  ImageDimension image_width_1 = 6000;
  ImageDimension image_height_1 = 4000;
  Scalar pixel_size_1 = 0.0000039;
  size_t number_of_points_1 = 2000;
  Scalar lateral_overlap_ratio_1 = 0.6;
  Scalar longitudinal_overlap_ratio_1 = 0.8;
  Scalar scene_max_height_1 = 50;
  Scalar camera_height_stddev_1 = 2;
  Scalar camera_planar_stddev_1 = 2;
  Scalar camera_rotation_stddev_1 = 10;
  FlightGenerator flight_generator_1(
    focal_length_in_metre_1,
    number_of_strips_1,
    number_of_cameras_in_strip_1,
    ground_resolution_1,
    image_width_1,
    image_height_1,
    pixel_size_1,
    number_of_points_1,
    lateral_overlap_ratio_1,
    longitudinal_overlap_ratio_1,
    scene_max_height_1,
    camera_height_stddev_1,
    camera_planar_stddev_1,
    camera_rotation_stddev_1);
  flight_generators.push_back(flight_generator_1);

  Scalar focal_length_in_metre_2 = 0.019;
  size_t number_of_strips_2 = 10;
  size_t number_of_cameras_in_strip_2 = 30;
  Scalar ground_resolution_2 = 0.1;
  ImageDimension image_width_2 = 6000;
  ImageDimension image_height_2 = 4000;
  Scalar pixel_size_2 = 0.0000039;
  size_t number_of_points_2 = 2000;
  Scalar lateral_overlap_ratio_2 = 0.6;
  Scalar longitudinal_overlap_ratio_2 = 0.8;
  Scalar scene_max_height_2 = 50;
  Scalar camera_height_stddev_2 = 2;
  Scalar camera_planar_stddev_2 = 2;
  Scalar camera_rotation_stddev_2 = 10;
  FlightGenerator flight_generator_2(
    focal_length_in_metre_2,
    number_of_strips_2,
    number_of_cameras_in_strip_2,
    ground_resolution_2,
    image_width_2,
    image_height_2,
    pixel_size_2,
    number_of_points_2,
    lateral_overlap_ratio_2,
    longitudinal_overlap_ratio_2,
    scene_max_height_2,
    camera_height_stddev_2,
    camera_planar_stddev_2,
    camera_rotation_stddev_2);
  flight_generators.push_back(flight_generator_2);

  Test test(0.8, 0.2, 60, 10, 15, flight_generators, 5000, "simple_test");
  ASSERT_EQ(0, test.Test());
}

}
