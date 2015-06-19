#include <gtest/gtest.h>

#include "hs_sfm/sfm_file_io/scene_ply_saver.hpp"
#include "hs_sfm/synthetic/multiple_flight_generator.hpp"
#include "hs_sfm/sfm_pipeline/point_cloud_norm_calculator.hpp"

namespace
{

TEST(TestPointCloudNormCalculator, SyntheticSimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef hs::sfm::pipeline::PointCloudNormCalculator<Scalar> Calculator;
  typedef Calculator::PointContainer PointContainer;
  typedef Calculator::Point Point;
  typedef Calculator::IntrinsicParams IntrinsicParams;
  typedef Calculator::ExtrinsicParams ExtrinsicParams;
  typedef Calculator::CameraParams CameraParams;
  typedef Calculator::CameraParamsContainer CameraParamsContainer;

  typedef hs::sfm::synthetic::MultipleFlightGenerator<Scalar, ImageDimension>
          Generator;
  typedef Generator::FlightGenerator FlightGenerator;
  typedef Generator::FlightGeneratorContainer FlightGeneratorContainer;
  typedef Generator::Image Image;
  typedef Generator::ImageContainer ImageContainer;

  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef Generator::ExtrinsicParamsContainer ExtrinsicParamsContainer;

  typedef hs::sfm::fileio::ScenePLYSaver<Scalar, ImageDimension> Saver;

  FlightGeneratorContainer flight_generators;
  IntrinsicParamsContainer intrinsic_params_set_true;
  IntrinsicParamsContainer intrinsic_params_set_initial;

  Scalar focal_length_in_metre_0 = 0.018858358970276164;
  size_t number_of_strips_0 = 3;
  size_t number_of_cameras_in_strip_0 = 30;
  Scalar ground_resolution_0 = 0.1;
  ImageDimension image_width_0 = 6000;
  ImageDimension image_height_0 = 4000;
  Scalar pixel_size_0 = 0.0000039;
  size_t number_of_points_0 = 100000;
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
                                     3000-42.4095312016,
                                     2000+31.699212823,
                                     1,
                                     -0.0050490462006048831,
                                     0.031293804298609881,
                                     -0.030794820960442223,
                                     -0.00055376548320189127,
                                     -0.00049877717768381476);
  intrinsic_params_set_true.push_back(intrinsic_params_0);

  Scalar focal_length_in_metre_1 = 0.02995452167701055;
  size_t number_of_strips_1 = 3;
  size_t number_of_cameras_in_strip_1 = 30;
  Scalar ground_resolution_1 = 0.1;
  ImageDimension image_width_1 = 6000;
  ImageDimension image_height_1 = 4000;
  Scalar pixel_size_1 = 0.0000039;
  size_t number_of_points_1 = 100000;
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
                                     3000-21.669436058,
                                     2000+44.8644764322,
                                     1,
                                     -0.02529179096221609,
                                     0.23762413973445157,
                                     -0.64208397668697237,
                                     -0.0020605099808780948,
                                     -0.00028706423764766859);
  intrinsic_params_set_true.push_back(intrinsic_params_1);

  Scalar focal_length_in_metre_2 = 0.019056097774998712;
  size_t number_of_strips_2 = 3;
  size_t number_of_cameras_in_strip_2 = 30;
  Scalar ground_resolution_2 = 0.1;
  ImageDimension image_width_2 = 6000;
  ImageDimension image_height_2 = 4000;
  Scalar pixel_size_2 = 0.0000039;
  size_t number_of_points_2 = 100000;
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
                                     3000-35.2052431556,
                                     2000+16.4262220759,
                                     1,
                                     -0.10316088386868619,
                                     0.13520490482776426,
                                     -0.05489235547426094,
                                     4.1434720317373253e-006,
                                     -0.00025018439997095336);
  intrinsic_params_set_true.push_back(intrinsic_params_2);

  Scalar flight_longitudinal_overlap_ratio = 0.9;
  Scalar flight_lateral_overlap_ratio = 0.2;
  Scalar north_west_angle = 0.0;
  Scalar north_west_angle_stddev = 5.0;
  Scalar offset_stddev = 15.0;

  Generator generator(flight_longitudinal_overlap_ratio,
                      flight_lateral_overlap_ratio,
                      north_west_angle,
                      north_west_angle_stddev,
                      offset_stddev,
                      flight_generators);

  size_t number_of_points = 600000;
  PointContainer points;
  ASSERT_EQ(0, generator.GeneratePoints(number_of_points, points));

  CameraParamsContainer camera_params_set;
  IntrinsicParamsContainer intrinsic_params_set_out;
  ExtrinsicParamsContainer extrinsic_params_set_out;
  ImageContainer images_out;
  for (size_t i = 0; i < flight_generators.size(); i++)
  {
    ExtrinsicParamsContainer extrinsic_params_set_flight;
    ImageContainer images_flight;
    ASSERT_EQ(0, generator.GenerateExtrinsicParamsContainer(
                   i, extrinsic_params_set_flight, images_flight));
    for (size_t j = 0; j < extrinsic_params_set_flight.size(); j++)
    {
      CameraParams camera;
      camera.intrinsic_params = intrinsic_params_set_true[i];
      camera.extrinsic_params = extrinsic_params_set_flight[j];
      camera.image_width = images_flight[j].m_width;
      camera.image_height = images_flight[j].m_height;
      camera_params_set.push_back(camera);
      intrinsic_params_set_out.push_back(intrinsic_params_set_true[i]);
      extrinsic_params_set_out.push_back(extrinsic_params_set_flight[j]);
      images_out.push_back(images_flight[j]);
    }
  }

  Saver saver(Scalar(5));

  Calculator calculator;
  PointContainer norms;
  ASSERT_EQ(0, calculator(camera_params_set, points, norms));

  ASSERT_EQ(0, saver(std::string("point_cloud_norms_with_cameras.ply"),
                     intrinsic_params_set_out, extrinsic_params_set_out,
                     images_out, points, &norms));

  Point up_vector;
  up_vector << Scalar(0),
               Scalar(0),
               Scalar(1);
  ASSERT_EQ(0, calculator(points, up_vector, norms));
  ASSERT_EQ(0, saver(std::string("point_cloud_norms_with_up_vector.ply"),
                     intrinsic_params_set_out, extrinsic_params_set_out,
                     images_out, points, &norms));

}

}