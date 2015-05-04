#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "hs_sfm/sfm_file_io/keyset_saver.hpp"
#include "hs_sfm/sfm_file_io/tracks_saver.hpp"
#include "hs_sfm/sfm_file_io/intrinsic_params_set_saver.hpp"
#include "hs_sfm/sfm_file_io/extrinsic_params_set_saver.hpp"
#include "hs_sfm/sfm_file_io/points_saver.hpp"
#include "hs_sfm/sfm_file_io/object_index_saver.hpp"
#include "hs_sfm/incremental/scene_expandor.hpp"

#include "synthetic_data_generator.hpp"
#include "data_tester.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestSceneExpandor
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
  typedef hs::sfm::incremental::SyntheticDataGenerator<Scalar, ImageDimension>
          DataGenerator;

  typedef typename DataGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename DataGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename DataGenerator::ImageContainer ImageContainer;
  typedef typename DataGenerator::Point3DContainer Point3DContainer;
  typedef typename DataGenerator::KeysetContainer KeysetContainer;
  typedef typename DataGenerator::RMatrix RMatrix;
  typedef typename DataGenerator::Translate Translate;

  typedef hs::sfm::incremental::SceneExpandor<Scalar> Expandor;

  typedef hs::sfm::incremental::DataTester<Scalar> Tester;

public:
  TestSceneExpandor(
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
    Scalar camera_rotation_stddev,
    Scalar north_west_angle,
    Scalar outlier_ratio,
    Scalar key_stddev,
    const std::string& test_name,
    size_t number_of_gcps)
  : data_generator_(focal_length_in_metre,
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
                    outlier_ratio,
                    key_stddev),
    test_name_(test_name),
    number_of_gcps_(number_of_gcps) {}

  Err Test()
  {
    IntrinsicParamsContainer intrinsic_params_set;
    ExtrinsicParamsContainer extrinsic_params_set_absolute;
    ImageContainer images;
    Point3DContainer points_absolute;
    KeysetContainer keysets_true;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;

    if (data_generator_.GenerateAbsoluteScene(intrinsic_params_set,
                                              extrinsic_params_set_absolute,
                                              images,
                                              points_absolute,
                                              keysets_true,
                                              tracks,
                                              camera_views) != 0)
    {
      return -1;
    }

    hs::sfm::fileio::ExtrinsicParamsSetSaver<Scalar> extrinsic_set_saver;
    hs::sfm::fileio::PointsSaver<Scalar> points_saver;
    hs::sfm::fileio::TracksSaver tracks_saver;
    hs::sfm::fileio::IntrinsicParamsSetSaver<Scalar> intrinsic_set_saver;
    hs::sfm::fileio::ObjectIndexSaver object_index_saver;
    hs::sfm::fileio::KeysetSaver<Scalar> keyset_saver;

    std::string extrinsic_set_absolute_true_path =
      test_name_ + "_extrinsic_params_set_absolute_true.txt";
    if (extrinsic_set_saver(extrinsic_set_absolute_true_path,
                            extrinsic_params_set_absolute) != 0)
      return -1;
    std::string points_absolute_true_path =
      test_name_ + "_points_absolute_true.txt";
    if (points_saver(points_absolute_true_path,
                     points_absolute) != 0)
      return -1;

    ExtrinsicParamsContainer extrinsic_params_set_relative;
    Point3DContainer points_relative;
    RMatrix rotation_similar;
    Translate translate_similar;
    Scalar scale_similar;
    size_t camera_id_identity;
    size_t camera_id_relative;
    if (data_generator_.GenerateRelativeScene(tracks,
                                              camera_views,
                                              extrinsic_params_set_absolute,
                                              points_absolute,
                                              rotation_similar,
                                              translate_similar,
                                              scale_similar,
                                              camera_id_identity,
                                              camera_id_relative,
                                              extrinsic_params_set_relative,
                                              points_relative) != 0)
    {
      return -1;
    }

    KeysetContainer keysets_noised;
    if (data_generator_.GenerateNoisedKeysets(keysets_true,
                                              keysets_noised) != 0)
    {
      return -1;
    }
    for (size_t i = 0; i < keysets_noised.size(); i++)
    {
      std::stringstream ss;
      ss<<test_name_<<"_keyset_"<<i<<".txt";
      std::string keyset_path;
      ss>>keyset_path;
      if (keyset_saver(keyset_path, keysets_noised[i]) != 0) return -1;
    }

    KeysetContainer keysets_gcp_true;
    Point3DContainer gcps;
    hs::sfm::TrackContainer tracks_gcp;
    if (data_generator_.GenerateGCPData(intrinsic_params_set,
                                        extrinsic_params_set_absolute,
                                        number_of_gcps_,
                                        keysets_gcp_true,
                                        gcps,
                                        tracks_gcp) != 0)
    {
      return -1;
    }
    hs::sfm::ViewInfoIndexer view_info_indexer_gcp;
    view_info_indexer_gcp.SetViewInfoByTracks(tracks_gcp);
    std::string tracks_gcp_path = test_name_ + "_tracks_gcp.txt";
    if (tracks_saver(tracks_gcp_path, tracks_gcp, view_info_indexer_gcp) != 0)
    {
      return -1;
    }

    std::string gcps_path = test_name_ + "_gcps.txt";
    if (points_saver(gcps_path, gcps) != 0)
    {
      return -1;
    }

    KeysetContainer keysets_gcp_noised;
    if (data_generator_.GenerateNoisedKeysets(keysets_gcp_true,
                                              keysets_gcp_noised) != 0)
    {
      return -1;
    }
    for (size_t i = 0; i < keysets_gcp_noised.size(); i++)
    {
      std::stringstream ss;
      ss<<test_name_<<"_keyset_gcp_"<<i<<".txt";
      std::string keyset_gcp_path;
      ss>>keyset_gcp_path;
      if (keyset_saver(keyset_gcp_path, keysets_gcp_noised[i]) != 0) return -1;
    }

    ExtrinsicParamsContainer extrinsic_params_set_relative_estimate;
    Point3DContainer points_relative_estimate;
    hs::sfm::ObjectIndexMap image_extrinsic_map;
    hs::sfm::ObjectIndexMap track_point_map;
    if (data_generator_.GenerateInitialScene(
          keysets_noised,
          intrinsic_params_set,
          tracks,
          camera_views,
          extrinsic_params_set_relative,
          points_relative,
          camera_id_identity,
          camera_id_relative,
          extrinsic_params_set_relative_estimate,
          points_relative_estimate,
          image_extrinsic_map,
          track_point_map) != 0)
    {
      return -1;
    }

    Expandor expandor(9, 4 * data_generator_.key_stddev());
    hs::sfm::ViewInfoIndexer view_info_indexer;
    if (expandor(keysets_noised,
                 intrinsic_params_set,
                 tracks,
                 extrinsic_params_set_relative_estimate,
                 image_extrinsic_map,
                 points_relative_estimate,
                 track_point_map,
                 view_info_indexer) != 0)
    {
      return -1;
    }

    Tester tester;
    if (tester.TestReprojectiveError(
          keysets_noised,
          intrinsic_params_set,
          tracks,
          image_extrinsic_map,
          track_point_map,
          view_info_indexer,
          extrinsic_params_set_relative_estimate,
          points_relative_estimate,
          data_generator_.key_stddev()) != 0)
    {
      return -1;
    }

    std::string tracks_path = test_name_ + "_tracks.txt";
    if (tracks_saver(tracks_path, tracks, view_info_indexer) != 0) return -1;

    std::string intrinsic_set_path =
      test_name_ + "_intrinsic_params_set.txt";
    if (intrinsic_set_saver(intrinsic_set_path,
                            intrinsic_params_set) != 0) return -1;

    std::string extrinsic_set_relative_estimate_path =
      test_name_ + "_extrinsic_params_set_relative_estimate.txt";
    if (extrinsic_set_saver(extrinsic_set_relative_estimate_path,
                            extrinsic_params_set_relative_estimate) != 0)
      return -1;

    std::string points_relative_estimate_path =
      test_name_ + "_points_relative_estimate.txt";
    if (points_saver(points_relative_estimate_path,
                     points_relative_estimate) != 0) return -1;

    std::string track_point_map_path = test_name_ + "_track_point_map.txt";
    if (object_index_saver(track_point_map_path,
                           track_point_map) != 0) return -1;

    std::string image_extrinsic_map_path =
      test_name_ + "_image_extrinsic_map.txt";
    if (object_index_saver(image_extrinsic_map_path,
                           image_extrinsic_map) != 0) return -1;

    return 0;
  }

private:
  DataGenerator data_generator_;
  std::string test_name_;
  size_t number_of_gcps_;
};

TEST(TestSceneExpandor, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestSceneExpandor<Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 3;
  size_t number_of_cameras_in_strips = 3;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 200;
  Scalar lateral_overlap_ratio = 0.7;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 100;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 5;
  Scalar north_west_angle = 60;
  Scalar outlier_ratio = 0.01;
  Scalar key_stddev = 2.0;
  std::string test_name = "small_data";
  size_t number_of_gcps = 10;

  Test test(focal_length_in_metre,
            number_of_strips,
            number_of_cameras_in_strips,
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
            outlier_ratio,
            key_stddev,
            test_name,
            number_of_gcps);

  ASSERT_EQ(0, test.Test());
}

TEST(TestSceneExpandor, BigDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestSceneExpandor<Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 7;
  size_t number_of_cameras_in_strips = 15;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 20000;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 100;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 5;
  Scalar north_west_angle = 60;
  Scalar outlier_ratio = 0.0;
  Scalar key_stddev = 1;
  std::string test_name = "big_data";
  size_t number_of_gcps = 10;

  Test test(focal_length_in_metre,
            number_of_strips,
            number_of_cameras_in_strips,
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
            outlier_ratio,
            key_stddev,
            test_name,
            number_of_gcps);

  ASSERT_EQ(0, test.Test());
}

}
