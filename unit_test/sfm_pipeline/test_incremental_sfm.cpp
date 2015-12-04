#include <gtest/gtest.h>

#include "data_tester.hpp"
#include "sfm_pipeline_tester.hpp"
#include "synthetic_data_generator.hpp"
#include "real_data_generator.hpp"
#include "real_synthetic_data_generator.hpp"

#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"
#include "hs_sfm/sfm_pipeline/incremental_sfm.hpp"
#include "hs_sfm/sfm_pipeline/gcp_similar_transform_estimator.hpp"
#include "hs_sfm/triangulate/multiple_view_maximum_likelihood_estimator.hpp"
#include "hs_sfm/sfm_pipeline/bundle_adjustment_gcp_constrained_optimizor.hpp"
#define DEBUG_TMP 1
#if DEBUG_TMP
#include "hs_sfm/sfm_utility/debug_tmp.hpp"
#include "hs_sfm/sfm_file_io/scene_ply_saver.hpp"
#include <fstream>
#include <iomanip>
#endif

namespace
{

TEST(TestIncrementalSFM, Synthetic120ImagesTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef hs::sfm::pipeline::IncrementalSFM<Scalar> SFMPipeline;	//主类
  typedef hs::sfm::pipeline::SFMPipelineTester<SFMPipeline> Tester;	//Tester容器
  typedef Tester::IntrinsicParams IntrinsicParams;	//注意此处直到下方的PointContainer所使用的是 IncrementalSFM<Scalar> 暴露的对应类型
  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Tester::ExtrinsicParams ExtrinsicParams;
  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef Tester::KeysetContainer KeysetContainer;
  typedef Tester::PointContainer PointContainer;
  typedef Tester::MatchContainer MatchContainer;
  typedef Tester::TrackContainer TrackContainer;

  typedef hs::sfm::pipeline::SyntheticDataGenerator<Scalar, ImageDimension>
          Generator;

  typedef Generator::FlightGenerator FlightGenerator;
  typedef Generator::FlightGeneratorContainer FlightGeneratorContainer;
  typedef Generator::Image Image;
  typedef Generator::ImageContainer ImageContainer;

  FlightGeneratorContainer flight_generators;
  IntrinsicParamsContainer intrinsic_params_set_true;
  IntrinsicParamsContainer intrinsic_params_set_initial;

  Scalar focal_length_in_metre_0 = 0.018858358970276164;
  size_t number_of_strips_0 = 2;
  size_t number_of_cameras_in_strip_0 = 30;
  Scalar ground_resolution_0 = 0.1;
  ImageDimension image_width_0 = 6000;
  ImageDimension image_height_0 = 4000;
  Scalar pixel_size_0 = 0.0000039;
  size_t number_of_points_0 = 10000;
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
  intrinsic_params_set_initial.push_back(
    IntrinsicParams(4871.79487179487179,
                    0,
                    Scalar(image_width_0 / 2),
                    Scalar(image_height_0 / 2)));

  Scalar focal_length_in_metre_1 = 0.02995452167701055;
  size_t number_of_strips_1 = 2;
  size_t number_of_cameras_in_strip_1 = 30;
  Scalar ground_resolution_1 = 0.1;
  ImageDimension image_width_1 = 6000;
  ImageDimension image_height_1 = 4000;
  Scalar pixel_size_1 = 0.0000039;
  size_t number_of_points_1 = 10000;
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
  //flight_generators.push_back(flight_generator_1);
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
  //intrinsic_params_set_true.push_back(intrinsic_params_1);
  //intrinsic_params_set_initial.push_back(
  //  IntrinsicParams(7692.30769230769231,
  //                  0,
  //                  Scalar(image_width_1 / 2),
  //                  Scalar(image_height_1 / 2)));

  Scalar focal_length_in_metre_2 = 0.019056097774998712;
  size_t number_of_strips_2 = 2;
  size_t number_of_cameras_in_strip_2 = 30;
  Scalar ground_resolution_2 = 0.1;
  ImageDimension image_width_2 = 6000;
  ImageDimension image_height_2 = 4000;
  Scalar pixel_size_2 = 0.0000039;
  size_t number_of_points_2 = 10000;
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
  intrinsic_params_set_initial.push_back(
    IntrinsicParams(4871.79487179487179,
                    0,
                    Scalar(image_width_2 / 2),
                    Scalar(image_height_2 / 2)));

  Scalar flight_longitudinal_overlap_ratio = 0.9;
  Scalar flight_lateral_overlap_ratio = 0.2;
  Scalar north_west_angle = 0.0;
  Scalar north_west_angle_stddev = 5.0;
  Scalar offset_stddev = 15.0;
  size_t number_of_points = 50000;
  size_t number_of_gcps = 10;
  size_t number_of_check_points = 500;
  //TODO:Outlier test needed!
  Scalar outlier_ratio = 0.1;
  Scalar key_stddev = 1.0;

  Generator generator(flight_longitudinal_overlap_ratio,
                      flight_lateral_overlap_ratio,
                      north_west_angle,
                      north_west_angle_stddev,
                      offset_stddev,
                      flight_generators,
                      outlier_ratio,
                      key_stddev,
                      intrinsic_params_set_true,
                      number_of_points,
                      number_of_gcps,
                      number_of_check_points);

  std::vector<size_t> image_intrinsic_map;
  hs::sfm::MatchContainer matches;
  KeysetContainer keysets_noised;
  PointContainer gcps;
  hs::sfm::TrackContainer tracks_gcp;
  KeysetContainer keysets_gcp_noised;
  PointContainer check_points;
  hs::sfm::TrackContainer tracks_check_point;
  KeysetContainer keysets_check_point_noised;
  ExtrinsicParamsContainer extrinsic_params_set_absolute_true;
  PointContainer points_absolute_true;
  ASSERT_EQ(0, generator.Generate(image_intrinsic_map,
                                  matches,
                                  keysets_noised,
                                  gcps,
                                  tracks_gcp,
                                  keysets_gcp_noised,
                                  check_points,
                                  tracks_check_point,
                                  keysets_check_point_noised,
                                  extrinsic_params_set_absolute_true,
                                  points_absolute_true));

  std::string test_name = "synthetic_data_120_images";

  SFMPipeline sfm_pipeline(100, 8, 2, 7);
  Tester tester(sfm_pipeline, 0.1, key_stddev, test_name);
  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
                      image_intrinsic_map,
                      matches,
                      keysets_noised,
                      gcps,
                      tracks_gcp,
                      keysets_gcp_noised,
                      check_points,
                      tracks_check_point,
                      keysets_check_point_noised,
                      intrinsic_params_set_true,
                      extrinsic_params_set_absolute_true,
                      points_absolute_true));

}

//TEST(TestIncrementalSFM, Synthetic240ImagesTest)
//{
//  typedef double Scalar;
//  typedef size_t ImageDimension;
//  typedef hs::sfm::pipeline::IncrementalSFM<Scalar> SFMPipeline;
//  typedef hs::sfm::pipeline::SFMPipelineTester<SFMPipeline> Tester;
//  typedef Tester::IntrinsicParams IntrinsicParams;
//  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
//  typedef Tester::ExtrinsicParams ExtrinsicParams;
//  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
//  typedef Tester::KeysetContainer KeysetContainer;
//  typedef Tester::PointContainer PointContainer;
//  typedef Tester::MatchContainer MatchContainer;
//  typedef Tester::TrackContainer TrackContainer;
//
//  typedef hs::sfm::pipeline::SyntheticDataGenerator<Scalar, ImageDimension>
//          Generator;
//
//  typedef Generator::FlightGenerator FlightGenerator;
//  typedef Generator::FlightGeneratorContainer FlightGeneratorContainer;
//  typedef Generator::Image Image;
//  typedef Generator::ImageContainer ImageContainer;
//
//  FlightGeneratorContainer flight_generators;
//  IntrinsicParamsContainer intrinsic_params_set_true;
//  IntrinsicParamsContainer intrinsic_params_set_initial;
//
//  Scalar focal_length_in_metre_0 = 0.018858358970276164;
//  size_t number_of_strips_0 = 2;
//  size_t number_of_cameras_in_strip_0 = 60;
//  Scalar ground_resolution_0 = 0.1;
//  ImageDimension image_width_0 = 6000;
//  ImageDimension image_height_0 = 4000;
//  Scalar pixel_size_0 = 0.0000039;
//  size_t number_of_points_0 = 10000;
//  Scalar lateral_overlap_ratio_0 = 0.6;
//  Scalar longitudinal_overlap_ratio_0 = 0.8;
//  Scalar scene_max_height_0 = 50;
//  Scalar camera_height_stddev_0 = 2;
//  Scalar camera_planar_stddev_0 = 2;
//  Scalar camera_rotation_stddev_0 = 10;
//  FlightGenerator flight_generator_0(
//    focal_length_in_metre_0,
//    number_of_strips_0,
//    number_of_cameras_in_strip_0,
//    ground_resolution_0,
//    image_width_0,
//    image_height_0,
//    pixel_size_0,
//    number_of_points_0,
//    lateral_overlap_ratio_0,
//    longitudinal_overlap_ratio_0,
//    scene_max_height_0,
//    camera_height_stddev_0,
//    camera_planar_stddev_0,
//    camera_rotation_stddev_0);
//  flight_generators.push_back(flight_generator_0);
//  IntrinsicParams intrinsic_params_0(focal_length_in_metre_0 / pixel_size_0,
//                                     0,
//                                     3000-42.4095312016,
//                                     2000+31.699212823,
//                                     1,
//                                     -0.0050490462006048831,
//                                     0.031293804298609881,
//                                     -0.030794820960442223,
//                                     -0.00055376548320189127,
//                                     -0.00049877717768381476);
//  intrinsic_params_set_true.push_back(intrinsic_params_0);
//  intrinsic_params_set_initial.push_back(
//    IntrinsicParams(4871.79487179487179,
//                    0,
//                    Scalar(image_width_0 / 2),
//                    Scalar(image_height_0 / 2)));
//
//  Scalar focal_length_in_metre_1 = 0.02995452167701055;
//  size_t number_of_strips_1 = 2;
//  size_t number_of_cameras_in_strip_1 = 60;
//  Scalar ground_resolution_1 = 0.1;
//  ImageDimension image_width_1 = 6000;
//  ImageDimension image_height_1 = 4000;
//  Scalar pixel_size_1 = 0.0000039;
//  size_t number_of_points_1 = 10000;
//  Scalar lateral_overlap_ratio_1 = 0.6;
//  Scalar longitudinal_overlap_ratio_1 = 0.8;
//  Scalar scene_max_height_1 = 50;
//  Scalar camera_height_stddev_1 = 2;
//  Scalar camera_planar_stddev_1 = 2;
//  Scalar camera_rotation_stddev_1 = 10;
//  FlightGenerator flight_generator_1(
//    focal_length_in_metre_1,
//    number_of_strips_1,
//    number_of_cameras_in_strip_1,
//    ground_resolution_1,
//    image_width_1,
//    image_height_1,
//    pixel_size_1,
//    number_of_points_1,
//    lateral_overlap_ratio_1,
//    longitudinal_overlap_ratio_1,
//    scene_max_height_1,
//    camera_height_stddev_1,
//    camera_planar_stddev_1,
//    camera_rotation_stddev_1);
//  //flight_generators.push_back(flight_generator_1);
//  IntrinsicParams intrinsic_params_1(focal_length_in_metre_1 / pixel_size_1,
//                                     0,
//                                     3000-21.669436058,
//                                     2000+44.8644764322,
//                                     1,
//                                     -0.02529179096221609,
//                                     0.23762413973445157,
//                                     -0.64208397668697237,
//                                     -0.0020605099808780948,
//                                     -0.00028706423764766859);
//  //intrinsic_params_set_true.push_back(intrinsic_params_1);
//  //intrinsic_params_set_initial.push_back(
//  //  IntrinsicParams(7692.30769230769231,
//  //                  0,
//  //                  Scalar(image_width_1 / 2),
//  //                  Scalar(image_height_1 / 2)));
//
//  Scalar focal_length_in_metre_2 = 0.019056097774998712;
//  size_t number_of_strips_2 = 2;
//  size_t number_of_cameras_in_strip_2 = 60;
//  Scalar ground_resolution_2 = 0.1;
//  ImageDimension image_width_2 = 6000;
//  ImageDimension image_height_2 = 4000;
//  Scalar pixel_size_2 = 0.0000039;
//  size_t number_of_points_2 = 10000;
//  Scalar lateral_overlap_ratio_2 = 0.6;
//  Scalar longitudinal_overlap_ratio_2 = 0.8;
//  Scalar scene_max_height_2 = 50;
//  Scalar camera_height_stddev_2 = 2;
//  Scalar camera_planar_stddev_2 = 2;
//  Scalar camera_rotation_stddev_2 = 10;
//  FlightGenerator flight_generator_2(
//    focal_length_in_metre_2,
//    number_of_strips_2,
//    number_of_cameras_in_strip_2,
//    ground_resolution_2,
//    image_width_2,
//    image_height_2,
//    pixel_size_2,
//    number_of_points_2,
//    lateral_overlap_ratio_2,
//    longitudinal_overlap_ratio_2,
//    scene_max_height_2,
//    camera_height_stddev_2,
//    camera_planar_stddev_2,
//    camera_rotation_stddev_2);
//  flight_generators.push_back(flight_generator_2);
//  IntrinsicParams intrinsic_params_2(focal_length_in_metre_2 / pixel_size_2,
//                                     0,
//                                     3000-35.2052431556,
//                                     2000+16.4262220759,
//                                     1,
//                                     -0.10316088386868619,
//                                     0.13520490482776426,
//                                     -0.05489235547426094,
//                                     4.1434720317373253e-006,
//                                     -0.00025018439997095336);
//  intrinsic_params_set_true.push_back(intrinsic_params_2);
//  intrinsic_params_set_initial.push_back(
//    IntrinsicParams(4871.79487179487179,
//                    0,
//                    Scalar(image_width_2 / 2),
//                    Scalar(image_height_2 / 2)));
//
//  Scalar flight_longitudinal_overlap_ratio = 0.9;
//  Scalar flight_lateral_overlap_ratio = 0.2;
//  Scalar north_west_angle = 0.0;
//  Scalar north_west_angle_stddev = 5.0;
//  Scalar offset_stddev = 15.0;
//  size_t number_of_points = 100000;
//  size_t number_of_gcps = 10;
//  size_t number_of_check_points = 500;
//  //TODO:Outlier test needed!
//  Scalar outlier_ratio = 0.1;
//  Scalar key_stddev = 1.0;
//
//  Generator generator(flight_longitudinal_overlap_ratio,
//                      flight_lateral_overlap_ratio,
//                      north_west_angle,
//                      north_west_angle_stddev,
//                      offset_stddev,
//                      flight_generators,
//                      outlier_ratio,
//                      key_stddev,
//                      intrinsic_params_set_true,
//                      number_of_points,
//                      number_of_gcps,
//                      number_of_check_points);
//
//  std::vector<size_t> image_intrinsic_map;
//  hs::sfm::MatchContainer matches;
//  KeysetContainer keysets_noised;
//  PointContainer gcps;
//  hs::sfm::TrackContainer tracks_gcp;
//  KeysetContainer keysets_gcp_noised;
//  PointContainer check_points;
//  hs::sfm::TrackContainer tracks_check_point;
//  KeysetContainer keysets_check_point_noised;
//  ExtrinsicParamsContainer extrinsic_params_set_absolute_true;
//  PointContainer points_absolute_true;
//  ASSERT_EQ(0, generator.Generate(image_intrinsic_map,
//                                  matches,
//                                  keysets_noised,
//                                  gcps,
//                                  tracks_gcp,
//                                  keysets_gcp_noised,
//                                  check_points,
//                                  tracks_check_point,
//                                  keysets_check_point_noised,
//                                  extrinsic_params_set_absolute_true,
//                                  points_absolute_true));
//
//  std::string test_name = "synthetic_data_240_images";
//
//  SFMPipeline sfm_pipeline(100, 8, 2, 7);
//  Tester tester(sfm_pipeline, 0.1, key_stddev, test_name);
//  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
//                      image_intrinsic_map,
//                      matches,
//                      keysets_noised,
//                      gcps,
//                      tracks_gcp,
//                      keysets_gcp_noised,
//                      check_points,
//                      tracks_check_point,
//                      keysets_check_point_noised,
//                      intrinsic_params_set_true,
//                      extrinsic_params_set_absolute_true,
//                      points_absolute_true));
//
//}

//TEST(TestIncrementalSFM, Synthetic480ImagesTest)
//{
//  typedef double Scalar;
//  typedef size_t ImageDimension;
//  typedef hs::sfm::pipeline::IncrementalSFM<Scalar> SFMPipeline;
//  typedef hs::sfm::pipeline::SFMPipelineTester<SFMPipeline> Tester;
//  typedef Tester::IntrinsicParams IntrinsicParams;
//  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
//  typedef Tester::ExtrinsicParams ExtrinsicParams;
//  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
//  typedef Tester::KeysetContainer KeysetContainer;
//  typedef Tester::PointContainer PointContainer;
//  typedef Tester::MatchContainer MatchContainer;
//  typedef Tester::TrackContainer TrackContainer;
//
//  typedef hs::sfm::pipeline::SyntheticDataGenerator<Scalar, ImageDimension>
//          Generator;
//
//  typedef Generator::FlightGenerator FlightGenerator;
//  typedef Generator::FlightGeneratorContainer FlightGeneratorContainer;
//  typedef Generator::Image Image;
//  typedef Generator::ImageContainer ImageContainer;
//
//  FlightGeneratorContainer flight_generators;
//  IntrinsicParamsContainer intrinsic_params_set_true;
//  IntrinsicParamsContainer intrinsic_params_set_initial;
//
//  Scalar focal_length_in_metre_0 = 0.018858358970276164;
//  size_t number_of_strips_0 = 4;
//  size_t number_of_cameras_in_strip_0 = 60;
//  Scalar ground_resolution_0 = 0.1;
//  ImageDimension image_width_0 = 6000;
//  ImageDimension image_height_0 = 4000;
//  Scalar pixel_size_0 = 0.0000039;
//  size_t number_of_points_0 = 10000;
//  Scalar lateral_overlap_ratio_0 = 0.6;
//  Scalar longitudinal_overlap_ratio_0 = 0.8;
//  Scalar scene_max_height_0 = 100;
//  Scalar camera_height_stddev_0 = 20;
//  Scalar camera_planar_stddev_0 = 20;
//  Scalar camera_rotation_stddev_0 = 60;
//  FlightGenerator flight_generator_0(
//    focal_length_in_metre_0,
//    number_of_strips_0,
//    number_of_cameras_in_strip_0,
//    ground_resolution_0,
//    image_width_0,
//    image_height_0,
//    pixel_size_0,
//    number_of_points_0,
//    lateral_overlap_ratio_0,
//    longitudinal_overlap_ratio_0,
//    scene_max_height_0,
//    camera_height_stddev_0,
//    camera_planar_stddev_0,
//    camera_rotation_stddev_0);
//  flight_generators.push_back(flight_generator_0);
//  IntrinsicParams intrinsic_params_0(focal_length_in_metre_0 / pixel_size_0,
//                                     0,
//                                     3000-42.4095312016,
//                                     2000+31.699212823,
//                                     1,
//                                     -0.0050490462006048831,
//                                     0.031293804298609881,
//                                     -0.030794820960442223,
//                                     -0.00055376548320189127,
//                                     -0.00049877717768381476);
//  intrinsic_params_set_true.push_back(intrinsic_params_0);
//  intrinsic_params_set_initial.push_back(
//    IntrinsicParams(4871.79487179487179,
//                    0,
//                    Scalar(image_width_0 / 2),
//                    Scalar(image_height_0 / 2)));
//
//  Scalar focal_length_in_metre_1 = 0.02995452167701055;
//  size_t number_of_strips_1 = 4;
//  size_t number_of_cameras_in_strip_1 = 60;
//  Scalar ground_resolution_1 = 0.1;
//  ImageDimension image_width_1 = 6000;
//  ImageDimension image_height_1 = 4000;
//  Scalar pixel_size_1 = 0.0000039;
//  size_t number_of_points_1 = 10000;
//  Scalar lateral_overlap_ratio_1 = 0.6;
//  Scalar longitudinal_overlap_ratio_1 = 0.8;
//  Scalar scene_max_height_1 = 100;
//  Scalar camera_height_stddev_1 = 20;
//  Scalar camera_planar_stddev_1 = 20;
//  Scalar camera_rotation_stddev_1 = 60;
//  FlightGenerator flight_generator_1(
//    focal_length_in_metre_1,
//    number_of_strips_1,
//    number_of_cameras_in_strip_1,
//    ground_resolution_1,
//    image_width_1,
//    image_height_1,
//    pixel_size_1,
//    number_of_points_1,
//    lateral_overlap_ratio_1,
//    longitudinal_overlap_ratio_1,
//    scene_max_height_1,
//    camera_height_stddev_1,
//    camera_planar_stddev_1,
//    camera_rotation_stddev_1);
//  //flight_generators.push_back(flight_generator_1);
//  IntrinsicParams intrinsic_params_1(focal_length_in_metre_1 / pixel_size_1,
//                                     0,
//                                     3000-21.669436058,
//                                     2000+44.8644764322,
//                                     1,
//                                     -0.02529179096221609,
//                                     0.23762413973445157,
//                                     -0.64208397668697237,
//                                     -0.0020605099808780948,
//                                     -0.00028706423764766859);
//  //intrinsic_params_set_true.push_back(intrinsic_params_1);
//  //intrinsic_params_set_initial.push_back(
//  //  IntrinsicParams(7692.30769230769231,
//  //                  0,
//  //                  Scalar(image_width_1 / 2),
//  //                  Scalar(image_height_1 / 2)));
//
//  Scalar focal_length_in_metre_2 = 0.019056097774998712;
//  size_t number_of_strips_2 = 4;
//  size_t number_of_cameras_in_strip_2 = 60;
//  Scalar ground_resolution_2 = 0.1;
//  ImageDimension image_width_2 = 6000;
//  ImageDimension image_height_2 = 4000;
//  Scalar pixel_size_2 = 0.0000039;
//  size_t number_of_points_2 = 10000;
//  Scalar lateral_overlap_ratio_2 = 0.6;
//  Scalar longitudinal_overlap_ratio_2 = 0.8;
//  Scalar scene_max_height_2 = 100;
//  Scalar camera_height_stddev_2 = 20;
//  Scalar camera_planar_stddev_2 = 20;
//  Scalar camera_rotation_stddev_2 = 60;
//  FlightGenerator flight_generator_2(
//    focal_length_in_metre_2,
//    number_of_strips_2,
//    number_of_cameras_in_strip_2,
//    ground_resolution_2,
//    image_width_2,
//    image_height_2,
//    pixel_size_2,
//    number_of_points_2,
//    lateral_overlap_ratio_2,
//    longitudinal_overlap_ratio_2,
//    scene_max_height_2,
//    camera_height_stddev_2,
//    camera_planar_stddev_2,
//    camera_rotation_stddev_2);
//  flight_generators.push_back(flight_generator_2);
//  IntrinsicParams intrinsic_params_2(focal_length_in_metre_2 / pixel_size_2,
//                                     0,
//                                     3000-35.2052431556,
//                                     2000+16.4262220759,
//                                     1,
//                                     -0.10316088386868619,
//                                     0.13520490482776426,
//                                     -0.05489235547426094,
//                                     4.1434720317373253e-006,
//                                     -0.00025018439997095336);
//  intrinsic_params_set_true.push_back(intrinsic_params_2);
//  intrinsic_params_set_initial.push_back(
//    IntrinsicParams(4871.79487179487179,
//                    0,
//                    Scalar(image_width_2 / 2),
//                    Scalar(image_height_2 / 2)));
//
//  Scalar flight_longitudinal_overlap_ratio = 0.9;
//  Scalar flight_lateral_overlap_ratio = 0.2;
//  Scalar north_west_angle = 0.0;
//  Scalar north_west_angle_stddev = 5.0;
//  Scalar offset_stddev = 15.0;
//  size_t number_of_points = 200000;
//  size_t number_of_gcps = 10;
//  size_t number_of_check_points = 500;
//  //TODO:Outlier test needed!
//  Scalar outlier_ratio = 0.1;
//  Scalar key_stddev = 1.0;
//
//  Generator generator(flight_longitudinal_overlap_ratio,
//                      flight_lateral_overlap_ratio,
//                      north_west_angle,
//                      north_west_angle_stddev,
//                      offset_stddev,
//                      flight_generators,
//                      outlier_ratio,
//                      key_stddev,
//                      intrinsic_params_set_true,
//                      number_of_points,
//                      number_of_gcps,
//                      number_of_check_points);
//
//  std::vector<size_t> image_intrinsic_map;
//  hs::sfm::MatchContainer matches;
//  KeysetContainer keysets_noised;
//  PointContainer gcps;
//  hs::sfm::TrackContainer tracks_gcp;
//  KeysetContainer keysets_gcp_noised;
//  PointContainer check_points;
//  hs::sfm::TrackContainer tracks_check_point;
//  KeysetContainer keysets_check_point_noised;
//  ExtrinsicParamsContainer extrinsic_params_set_absolute_true;
//  PointContainer points_absolute_true;
//  ASSERT_EQ(0, generator.Generate(image_intrinsic_map,
//                                  matches,
//                                  keysets_noised,
//                                  gcps,
//                                  tracks_gcp,
//                                  keysets_gcp_noised,
//                                  check_points,
//                                  tracks_check_point,
//                                  keysets_check_point_noised,
//                                  extrinsic_params_set_absolute_true,
//                                  points_absolute_true));
//
//  std::string test_name = "synthetic_data_480_images";
//
//  SFMPipeline sfm_pipeline(100, 8, 2, 7);
//  Tester tester(sfm_pipeline, 0.1, key_stddev, test_name);
//  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
//                      image_intrinsic_map,
//                      matches,
//                      keysets_noised,
//                      gcps,
//                      tracks_gcp,
//                      keysets_gcp_noised,
//                      check_points,
//                      tracks_check_point,
//                      keysets_check_point_noised,
//                      intrinsic_params_set_true,
//                      extrinsic_params_set_absolute_true,
//                      points_absolute_true));
//
//}

//TEST(TestIncrementalSFM, Synthetic960ImagesTest)
//{
//  typedef double Scalar;
//  typedef size_t ImageDimension;
//  typedef hs::sfm::pipeline::IncrementalSFM<Scalar> SFMPipeline;
//  typedef hs::sfm::pipeline::SFMPipelineTester<SFMPipeline> Tester;
//  typedef Tester::IntrinsicParams IntrinsicParams;
//  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
//  typedef Tester::ExtrinsicParams ExtrinsicParams;
//  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
//  typedef Tester::KeysetContainer KeysetContainer;
//  typedef Tester::PointContainer PointContainer;
//  typedef Tester::MatchContainer MatchContainer;
//  typedef Tester::TrackContainer TrackContainer;
//
//  typedef hs::sfm::pipeline::SyntheticDataGenerator<Scalar, ImageDimension>
//          Generator;
//
//  typedef Generator::FlightGenerator FlightGenerator;
//  typedef Generator::FlightGeneratorContainer FlightGeneratorContainer;
//  typedef Generator::Image Image;
//  typedef Generator::ImageContainer ImageContainer;
//
//  FlightGeneratorContainer flight_generators;
//  IntrinsicParamsContainer intrinsic_params_set_true;
//  IntrinsicParamsContainer intrinsic_params_set_initial;
//
//  Scalar focal_length_in_metre_0 = 0.018858358970276164;
//  size_t number_of_strips_0 = 8;
//  size_t number_of_cameras_in_strip_0 = 60;
//  Scalar ground_resolution_0 = 0.1;
//  ImageDimension image_width_0 = 6000;
//  ImageDimension image_height_0 = 4000;
//  Scalar pixel_size_0 = 0.0000039;
//  size_t number_of_points_0 = 10000;
//  Scalar lateral_overlap_ratio_0 = 0.6;
//  Scalar longitudinal_overlap_ratio_0 = 0.8;
//  Scalar scene_max_height_0 = 50;
//  Scalar camera_height_stddev_0 = 2;
//  Scalar camera_planar_stddev_0 = 2;
//  Scalar camera_rotation_stddev_0 = 10;
//  FlightGenerator flight_generator_0(
//    focal_length_in_metre_0,
//    number_of_strips_0,
//    number_of_cameras_in_strip_0,
//    ground_resolution_0,
//    image_width_0,
//    image_height_0,
//    pixel_size_0,
//    number_of_points_0,
//    lateral_overlap_ratio_0,
//    longitudinal_overlap_ratio_0,
//    scene_max_height_0,
//    camera_height_stddev_0,
//    camera_planar_stddev_0,
//    camera_rotation_stddev_0);
//  flight_generators.push_back(flight_generator_0);
//  IntrinsicParams intrinsic_params_0(focal_length_in_metre_0 / pixel_size_0,
//                                     0,
//                                     3000-42.4095312016,
//                                     2000+31.699212823,
//                                     1,
//                                     -0.0050490462006048831,
//                                     0.031293804298609881,
//                                     -0.030794820960442223,
//                                     -0.00055376548320189127,
//                                     -0.00049877717768381476);
//  intrinsic_params_set_true.push_back(intrinsic_params_0);
//  intrinsic_params_set_initial.push_back(
//    IntrinsicParams(4871.79487179487179,
//                    0,
//                    Scalar(image_width_0 / 2),
//                    Scalar(image_height_0 / 2)));
//
//  Scalar focal_length_in_metre_1 = 0.02995452167701055;
//  size_t number_of_strips_1 = 8;
//  size_t number_of_cameras_in_strip_1 = 60;
//  Scalar ground_resolution_1 = 0.1;
//  ImageDimension image_width_1 = 6000;
//  ImageDimension image_height_1 = 4000;
//  Scalar pixel_size_1 = 0.0000039;
//  size_t number_of_points_1 = 10000;
//  Scalar lateral_overlap_ratio_1 = 0.6;
//  Scalar longitudinal_overlap_ratio_1 = 0.8;
//  Scalar scene_max_height_1 = 50;
//  Scalar camera_height_stddev_1 = 2;
//  Scalar camera_planar_stddev_1 = 2;
//  Scalar camera_rotation_stddev_1 = 10;
//  FlightGenerator flight_generator_1(
//    focal_length_in_metre_1,
//    number_of_strips_1,
//    number_of_cameras_in_strip_1,
//    ground_resolution_1,
//    image_width_1,
//    image_height_1,
//    pixel_size_1,
//    number_of_points_1,
//    lateral_overlap_ratio_1,
//    longitudinal_overlap_ratio_1,
//    scene_max_height_1,
//    camera_height_stddev_1,
//    camera_planar_stddev_1,
//    camera_rotation_stddev_1);
//  //flight_generators.push_back(flight_generator_1);
//  IntrinsicParams intrinsic_params_1(focal_length_in_metre_1 / pixel_size_1,
//                                     0,
//                                     3000-21.669436058,
//                                     2000+44.8644764322,
//                                     1,
//                                     -0.02529179096221609,
//                                     0.23762413973445157,
//                                     -0.64208397668697237,
//                                     -0.0020605099808780948,
//                                     -0.00028706423764766859);
//  //intrinsic_params_set_true.push_back(intrinsic_params_1);
//  //intrinsic_params_set_initial.push_back(
//  //  IntrinsicParams(7692.30769230769231,
//  //                  0,
//  //                  Scalar(image_width_1 / 2),
//  //                  Scalar(image_height_1 / 2)));
//
//  Scalar focal_length_in_metre_2 = 0.019056097774998712;
//  size_t number_of_strips_2 = 8;
//  size_t number_of_cameras_in_strip_2 = 60;
//  Scalar ground_resolution_2 = 0.1;
//  ImageDimension image_width_2 = 6000;
//  ImageDimension image_height_2 = 4000;
//  Scalar pixel_size_2 = 0.0000039;
//  size_t number_of_points_2 = 10000;
//  Scalar lateral_overlap_ratio_2 = 0.6;
//  Scalar longitudinal_overlap_ratio_2 = 0.8;
//  Scalar scene_max_height_2 = 50;
//  Scalar camera_height_stddev_2 = 2;
//  Scalar camera_planar_stddev_2 = 2;
//  Scalar camera_rotation_stddev_2 = 10;
//  FlightGenerator flight_generator_2(
//    focal_length_in_metre_2,
//    number_of_strips_2,
//    number_of_cameras_in_strip_2,
//    ground_resolution_2,
//    image_width_2,
//    image_height_2,
//    pixel_size_2,
//    number_of_points_2,
//    lateral_overlap_ratio_2,
//    longitudinal_overlap_ratio_2,
//    scene_max_height_2,
//    camera_height_stddev_2,
//    camera_planar_stddev_2,
//    camera_rotation_stddev_2);
//  flight_generators.push_back(flight_generator_2);
//  IntrinsicParams intrinsic_params_2(focal_length_in_metre_2 / pixel_size_2,
//                                     0,
//                                     3000-35.2052431556,
//                                     2000+16.4262220759,
//                                     1,
//                                     -0.10316088386868619,
//                                     0.13520490482776426,
//                                     -0.05489235547426094,
//                                     4.1434720317373253e-006,
//                                     -0.00025018439997095336);
//  intrinsic_params_set_true.push_back(intrinsic_params_2);
//  intrinsic_params_set_initial.push_back(
//    IntrinsicParams(4871.79487179487179,
//                    0,
//                    Scalar(image_width_2 / 2),
//                    Scalar(image_height_2 / 2)));
//
//  Scalar flight_longitudinal_overlap_ratio = 0.9;
//  Scalar flight_lateral_overlap_ratio = 0.2;
//  Scalar north_west_angle = 0.0;
//  Scalar north_west_angle_stddev = 5.0;
//  Scalar offset_stddev = 15.0;
//  size_t number_of_points = 400000;
//  size_t number_of_gcps = 10;
//  size_t number_of_check_points = 500;
//  //TODO:Outlier test needed!
//  Scalar outlier_ratio = 0.2;
//  Scalar key_stddev = 1.0;
//
//  Generator generator(flight_longitudinal_overlap_ratio,
//                      flight_lateral_overlap_ratio,
//                      north_west_angle,
//                      north_west_angle_stddev,
//                      offset_stddev,
//                      flight_generators,
//                      outlier_ratio,
//                      key_stddev,
//                      intrinsic_params_set_true,
//                      number_of_points,
//                      number_of_gcps,
//                      number_of_check_points);
//
//  std::vector<size_t> image_intrinsic_map;
//  hs::sfm::MatchContainer matches;
//  KeysetContainer keysets_noised;
//  PointContainer gcps;
//  hs::sfm::TrackContainer tracks_gcp;
//  KeysetContainer keysets_gcp_noised;
//  PointContainer check_points;
//  hs::sfm::TrackContainer tracks_check_point;
//  KeysetContainer keysets_check_point_noised;
//  ExtrinsicParamsContainer extrinsic_params_set_absolute_true;
//  PointContainer points_absolute_true;
//  ASSERT_EQ(0, generator.Generate(image_intrinsic_map,
//                                  matches,
//                                  keysets_noised,
//                                  gcps,
//                                  tracks_gcp,
//                                  keysets_gcp_noised,
//                                  check_points,
//                                  tracks_check_point,
//                                  keysets_check_point_noised,
//                                  extrinsic_params_set_absolute_true,
//                                  points_absolute_true));
//
//  std::string test_name = "synthetic_data_960_images";
//
//  SFMPipeline sfm_pipeline(100, 8, 2, 7);
//  Tester tester(sfm_pipeline, 0.1, key_stddev, test_name);
//  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
//                      image_intrinsic_map,
//                      matches,
//                      keysets_noised,
//                      gcps,
//                      tracks_gcp,
//                      keysets_gcp_noised,
//                      check_points,
//                      tracks_check_point,
//                      keysets_check_point_noised,
//                      intrinsic_params_set_true,
//                      extrinsic_params_set_absolute_true,
//                      points_absolute_true));
//
//}

TEST(TestIncrementalSFM, Real242ImagesTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef hs::sfm::pipeline::IncrementalSFM<Scalar> SFMPipeline;
  typedef hs::sfm::pipeline::SFMPipelineTester<SFMPipeline> Tester;
  typedef Tester::IntrinsicParams IntrinsicParams;
  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Tester::KeysetContainer KeysetContainer;
  typedef Tester::PointContainer PointContainer;

  typedef hs::sfm::pipeline::RealDataGenerator<Scalar> Generator;

  typedef hs::sfm::MatchContainer MatchContainer;
  typedef hs::sfm::TrackContainer TrackContainer;

  std::string out_path =
    "../../test_data/sfm_pipeline/real_data_242_images/bundler.out";
  std::string gcp_path =
    "../../test_data/sfm_pipeline/real_data_242_images/gcp.xml";
  IntrinsicParamsContainer intrinsic_params_set_initial;
  intrinsic_params_set_initial.push_back(IntrinsicParams(4666.67,
                                                         0,
                                                         3000,
                                                         2000));
  std::vector<size_t> image_intrinsic_map(242);
  for (size_t i = 0; i < 242; i++)
  {
    image_intrinsic_map[i] = 0;
  }

  KeysetContainer keysets;
  MatchContainer matches;
  PointContainer gcps;
  TrackContainer tracks_gcp;
  KeysetContainer keysets_gcp;
  PointContainer check_points;
  TrackContainer tracks_check_point;
  KeysetContainer keysets_check_point;
  ASSERT_EQ(0, Generator::LoadGCPs(
    gcp_path,
    gcps, tracks_gcp, keysets_gcp,
    check_points, tracks_check_point, keysets_check_point));
  ASSERT_EQ(0, Generator::LoadBundlerOutFile(out_path, 6000, 4000,
                                             keysets, matches, false));

  Scalar key_stddev = Scalar(1);
  std::string test_name = "real_data_242_images";
  SFMPipeline sfm_pipeline(100, 8, 2, 7);
  Tester tester(sfm_pipeline, 0.1, key_stddev, test_name);
  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
                      image_intrinsic_map,
                      matches,
                      keysets,
                      gcps,
                      tracks_gcp,
                      keysets_gcp,
                      check_points,
                      tracks_check_point,
                      keysets_check_point));
}

TEST(TestIncrementalSFM, Real1067ImagesTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef hs::sfm::pipeline::IncrementalSFM<Scalar> SFMPipeline;
  typedef hs::sfm::pipeline::SFMPipelineTester<SFMPipeline> Tester;
  typedef Tester::IntrinsicParams IntrinsicParams;
  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Tester::KeysetContainer KeysetContainer;
  typedef Tester::PointContainer PointContainer;

  typedef hs::sfm::pipeline::RealDataGenerator<Scalar> Generator;

  typedef hs::sfm::MatchContainer MatchContainer;
  typedef hs::sfm::TrackContainer TrackContainer;

  std::string out_path =
    "../../sfm_pipeline/real_data_1067_images/bundler.out";
  std::string gcp_path =
    "../../sfm_pipeline/real_data_1067_images/gcp.xml";
  IntrinsicParamsContainer intrinsic_params_set_initial;
  intrinsic_params_set_initial.push_back(IntrinsicParams(7500.0,
                                                         0,
                                                         3000,
                                                         2000));
  intrinsic_params_set_initial.push_back(IntrinsicParams(7500.0,
                                                         0,
                                                         3000,
                                                         2000));
  std::vector<size_t> image_intrinsic_map(1067);
  for (size_t i = 0; i < 491; i++)
  {
    image_intrinsic_map[i] = 0;
  }
  for (size_t i = 491; i < 1067; i++)
  {
    image_intrinsic_map[i] = 1;
  }

  KeysetContainer keysets;
  MatchContainer matches;
  PointContainer gcps;
  TrackContainer tracks_gcp;
  KeysetContainer keysets_gcp;
  PointContainer check_points;
  TrackContainer tracks_check_point;
  KeysetContainer keysets_check_point;
  ASSERT_EQ(0, Generator::LoadGCPs(
    gcp_path,
    gcps, tracks_gcp, keysets_gcp,
    check_points, tracks_check_point, keysets_check_point));
  ASSERT_EQ(0, Generator::LoadBundlerOutFile(out_path, 6000, 4000,
                                             keysets, matches));

  Scalar key_stddev = Scalar(1);
  std::string test_name = "real_data_1067_images";
  SFMPipeline sfm_pipeline(100, 8, 2, 7);
  Tester tester(sfm_pipeline, 0.1, key_stddev, test_name);
  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
                      image_intrinsic_map,
                      matches,
                      keysets,
                      gcps,
                      tracks_gcp,
                      keysets_gcp,
                      check_points,
                      tracks_check_point,
                      keysets_check_point));
}

TEST(TestIncrementalSFM, RealSynthetic242ImagesTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef hs::sfm::pipeline::IncrementalSFM<Scalar> SFMPipeline;
  typedef hs::sfm::pipeline::SFMPipelineTester<SFMPipeline> Tester;
  typedef Tester::IntrinsicParams IntrinsicParams;
  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Tester::ExtrinsicParams ExtrinsicParams;
  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef Tester::KeysetContainer KeysetContainer;
  typedef Tester::PointContainer PointContainer;

  typedef hs::sfm::pipeline::RealSyntheticDataGenerator<Scalar> Generator;

  typedef hs::sfm::MatchContainer MatchContainer;
  typedef hs::sfm::TrackContainer TrackContainer;

  std::string out_path =
    "../../test_data/sfm_pipeline/real_data_242_images/bundler.out";
  std::string gcp_path =
    "../../test_data/sfm_pipeline/real_data_242_images/gcp.xml";
  IntrinsicParamsContainer intrinsic_params_set_initial;
  intrinsic_params_set_initial.push_back(IntrinsicParams(4666.67,
                                                         0,
                                                         3000,
                                                         2000));
  IntrinsicParamsContainer intrinsic_params_set_true;
  intrinsic_params_set_true.push_back(IntrinsicParams(4880.22,
                                                      0,
                                                      3026.11,
                                                      1992.36,
                                                      1,
                                                      -0.100423,
                                                      0.128487,
                                                      -0.0482081,
                                                      -3.17902e-05,
                                                      7.75629e-05));
  //intrinsic_params_set_true.push_back(IntrinsicParams(4666.67,
  //                                                    0,
  //                                                    3000,
  //                                                    2000));
  std::vector<size_t> image_intrinsic_map(242);
  for (size_t i = 0; i < 242; i++)
  {
    image_intrinsic_map[i] = 0;
  }

  KeysetContainer keysets;
  TrackContainer tracks;
  MatchContainer matches;
  PointContainer gcps;
  TrackContainer tracks_gcp;
  KeysetContainer keysets_gcp;
  PointContainer check_points;
  TrackContainer tracks_check_point;
  KeysetContainer keysets_check_point;
  PointContainer points_true;
  ExtrinsicParamsContainer extrinsic_params_set_true;
  Scalar key_stddev = Scalar(1);
  Scalar outlier_ratio = Scalar(0.1);

  Generator generator(outlier_ratio, key_stddev);
  generator.Generate(out_path, gcp_path,
                     intrinsic_params_set_true, image_intrinsic_map,
                     6000, 4000, tracks, matches,
                     points_true, extrinsic_params_set_true, keysets,
                     gcps, tracks_gcp, keysets_gcp,
                     check_points, tracks_check_point, keysets_check_point);

  std::string test_name = "real_synthetic_data_242_images";
  SFMPipeline sfm_pipeline(100, 8, 2, 7);
  Tester tester(sfm_pipeline, 0.1, key_stddev, test_name);
  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
                      image_intrinsic_map,
                      matches,
                      keysets,
                      gcps,
                      tracks_gcp,
                      keysets_gcp,
                      check_points,
                      tracks_check_point,
                      keysets_check_point,
                      intrinsic_params_set_true,
                      extrinsic_params_set_true,
                      points_true));
}

//TEST(TestIncrementalSFM, RealSynthetic1067ImagesTest)
//{
//  typedef double Scalar;
//  typedef size_t ImageDimension;
//  typedef hs::sfm::pipeline::IncrementalSFM<Scalar> SFMPipeline;
//  typedef hs::sfm::pipeline::SFMPipelineTester<SFMPipeline> Tester;
//  typedef Tester::IntrinsicParams IntrinsicParams;
//  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
//  typedef Tester::ExtrinsicParams ExtrinsicParams;
//  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
//  typedef Tester::KeysetContainer KeysetContainer;
//  typedef Tester::PointContainer PointContainer;
//
//  typedef hs::sfm::pipeline::RealSyntheticDataGenerator<Scalar> Generator;
//
//  typedef hs::sfm::MatchContainer MatchContainer;
//  typedef hs::sfm::TrackContainer TrackContainer;
//
//  std::string out_path =
//    "../../sfm_pipeline/real_data_1067_images/bundler.out";
//  std::string gcp_path =
//    "../../sfm_pipeline/real_data_1067_images/gcp.xml";
//  IntrinsicParamsContainer intrinsic_params_set_initial;
//  intrinsic_params_set_initial.push_back(IntrinsicParams(7500.00,
//                                                         0,
//                                                         3000,
//                                                         2000));
//  intrinsic_params_set_initial.push_back(IntrinsicParams(7500.00,
//                                                         0,
//                                                         3000,
//                                                         2000));
//  IntrinsicParamsContainer intrinsic_params_set_true;
//  intrinsic_params_set_true.push_back(IntrinsicParams(7680.93,
//                                                      0,
//                                                      3005.66,
//                                                      2002.71,
//                                                      1,
//                                                      -0.139239,
//                                                      0.264343,
//                                                      0.157404));
//  intrinsic_params_set_true.push_back(IntrinsicParams(7687.72,
//                                                      0,
//                                                      3001.87,
//                                                      2010.12,
//                                                      1,
//                                                      -0.138811,
//                                                      0.250986,
//                                                      0.20556,
//                                                      -0.00040312,
//                                                      0.000140978));
//
//  std::vector<size_t> image_intrinsic_map(1067);
//  for (size_t i = 0; i < 491; i++)
//  {
//    image_intrinsic_map[i] = 0;
//  }
//  for (size_t i = 491; i < 1067; i++)
//  {
//    image_intrinsic_map[i] = 1;
//  }
//
//  KeysetContainer keysets;
//  TrackContainer tracks;
//  MatchContainer matches;
//  PointContainer gcps;
//  TrackContainer tracks_gcp;
//  KeysetContainer keysets_gcp;
//  PointContainer check_points;
//  TrackContainer tracks_check_point;
//  KeysetContainer keysets_check_point;
//  PointContainer points_true;
//  ExtrinsicParamsContainer extrinsic_params_set_true;
//  Scalar key_stddev = Scalar(1);
//  Scalar outlier_ratio = Scalar(0.0);
//
//  Generator generator(outlier_ratio, key_stddev);
//  generator.Generate(out_path, gcp_path,
//                     intrinsic_params_set_true, image_intrinsic_map,
//                     6000, 4000, tracks, matches,
//                     points_true, extrinsic_params_set_true, keysets,
//                     gcps, tracks_gcp, keysets_gcp,
//                     check_points, tracks_check_point, keysets_check_point);
//
//  std::string test_name = "real_synthetic_data_1067_images";
//  SFMPipeline sfm_pipeline(100, 8, 2, 7);
//  Tester tester(sfm_pipeline, 0.1, key_stddev, test_name);
//  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
//                      image_intrinsic_map,
//                      matches,
//                      keysets,
//                      gcps,
//                      tracks_gcp,
//                      keysets_gcp,
//                      check_points,
//                      tracks_check_point,
//                      keysets_check_point,
//                      intrinsic_params_set_true,
//                      extrinsic_params_set_true,
//                      points_true));
//}

}
