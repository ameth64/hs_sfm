#include <iostream>

#include "gtest/gtest.h"

#include "hs_sfm/bundle_adjustment/camera_shared_synthetic_data_generator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"

namespace
{

TEST(TestCameraSharedVectorFunction, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef hs::sfm::ba::CameraSharedSyntheticDataGenerator<Scalar,
                                                          ImageDimension>
          DataGenerator;
  typedef DataGenerator::FlightGenerator FlightGenerator;
  typedef DataGenerator::FlightGeneratorContainer FlightGeneratorContainer;
  typedef DataGenerator::IntrinsicParams IntrinsicParams;
  typedef DataGenerator::IntrinsicParamsContainer IntrinsicParamsContainer;

  typedef hs::sfm::ba::CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef VectorFunction::Index Index;
  typedef VectorFunction::XVector XVector;
  typedef VectorFunction::YVector YVector;

  FlightGeneratorContainer flight_generators;
  IntrinsicParamsContainer intrinsic_params_set;

  Scalar focal_length_in_metre_0 = 0.018858358970276164;
  size_t number_of_strips_0 = 3;
  size_t number_of_cameras_in_strip_0 = 3;
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
  IntrinsicParams intrinsic_params_0(focal_length_in_metre_0 / pixel_size_0,
                                     0,
                                     -42.4095312016,
                                     -31.699212823,
                                     1,
                                     -0.0050490462006048831,
                                     0.031293804298609881,
                                     -0.030794820960442223,
                                     -0.00055376548320189127,
                                     -0.00049877717768381476);
  intrinsic_params_set.push_back(intrinsic_params_0);

  Scalar focal_length_in_metre_1 = 0.02995452167701055;
  size_t number_of_strips_1 = 4;
  size_t number_of_cameras_in_strip_1 = 3;
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
  IntrinsicParams intrinsic_params_1(focal_length_in_metre_1 / pixel_size_1,
                                     0,
                                     -21.669436058,
                                     -44.8644764322,
                                     1,
                                     -0.02529179096221609,
                                     0.23762413973445157,
                                     -0.64208397668697237,
                                     -0.0020605099808780948,
                                     -0.00028706423764766859);
  intrinsic_params_set.push_back(intrinsic_params_1);

  Scalar focal_length_in_metre_2 = 0.019056097774998712;
  size_t number_of_strips_2 = 2;
  size_t number_of_cameras_in_strip_2 = 5;
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
  IntrinsicParams intrinsic_params_2(focal_length_in_metre_2 / pixel_size_2,
                                     0,
                                     -35.2052431556,
                                     -16.4262220759,
                                     1,
                                     -0.10316088386868619,
                                     0.13520490482776426,
                                     -0.05489235547426094,
                                     4.1434720317373253e-006,
                                     -0.00025018439997095336);
  intrinsic_params_set.push_back(intrinsic_params_2);

  Scalar flight_longitudinal_overlap_ratio = 0.8;
  Scalar flight_lateral_overlap_ratio = 0.2;
  Scalar north_west_angle = 60.0;
  Scalar north_west_angle_stddev = 10.0;
  Scalar offset_stddev = 15.0;
  size_t number_of_points = 100;
  size_t number_of_planar_constrained_points = 4;
  size_t number_of_full_constrained_points = 6;
  size_t number_of_constrained_images = 10;
  size_t number_of_constrained_cameras = 2;

  DataGenerator data_generator(flight_longitudinal_overlap_ratio,
                               flight_lateral_overlap_ratio,
                               north_west_angle,
                               north_west_angle_stddev,
                               offset_stddev,
                               flight_generators,
                               number_of_points,
                               number_of_planar_constrained_points,
                               number_of_full_constrained_points,
                               number_of_constrained_images,
                               number_of_constrained_cameras,
                               intrinsic_params_set);

  VectorFunction vector_function;
  XVector x;
  YVector y_synthetic, y;
  Scalar threshold = 1e-8;

  ASSERT_EQ(0, data_generator(vector_function, x, y_synthetic));
  ASSERT_EQ(0, vector_function(x, y));
  ASSERT_EQ(true, y_synthetic.isApprox(y, threshold));
}

}
