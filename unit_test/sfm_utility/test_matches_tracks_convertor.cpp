#include <gtest/gtest.h>

#include <algorithm>

#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/synthetic/keyset_generator.hpp"
#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"

#include "hs_sfm/sfm_file_io/tracks_saver.hpp"

namespace
{

class TestMatchesTracksConvertor
{
public:
  typedef int Err;

  Err Test(hs::sfm::TrackContainer& tracks_true) const
  {
    hs::sfm::MatchesTracksConvertor convertor;
    hs::sfm::MatchContainer matches;
    if (convertor(tracks_true, matches) != 0) return -1;
    hs::sfm::TrackContainer tracks_estimate;
    if (convertor(matches, tracks_estimate) != 0) return -1;

    auto itr_track_true = tracks_true.begin();
    auto itr_track_true_end = tracks_true.end();
    for (; itr_track_true != itr_track_true_end; ++itr_track_true)
    {
      std::sort(itr_track_true->begin(), itr_track_true->end());
    }
    std::sort(tracks_true.begin(), tracks_true.end());

    auto itr_track_estimate = tracks_estimate.begin();
    auto itr_track_estimate_end = tracks_estimate.end();
    for (; itr_track_estimate != itr_track_estimate_end; ++itr_track_estimate)
    {
      std::sort(itr_track_estimate->begin(), itr_track_estimate->end());
    }
    std::sort(tracks_estimate.begin(), tracks_estimate.end());

    hs::sfm::fileio::TracksSaver saver;
    hs::sfm::ViewInfoIndexer view_info_indexer_true;;
    view_info_indexer_true.SetViewInfoByTracks(tracks_true);
    saver("tracks_true.txt", tracks_true, view_info_indexer_true);
    hs::sfm::ViewInfoIndexer view_info_indexer_estimate;
    view_info_indexer_estimate.SetViewInfoByTracks(tracks_estimate);
    saver("tracks_estimate.txt", tracks_estimate, view_info_indexer_estimate);

    if (tracks_true == tracks_estimate)
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }
};

TEST(TestMatchesTracksConvertor, SyntheticTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef hs::sfm::synthetic::SceneGenerator<Scalar, ImageDimension>
          SceneGenerator;
  typedef hs::sfm::synthetic::KeysetGenerator<Scalar, ImageDimension>
          KeysetGenerator;

  typedef SceneGenerator::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef SceneGenerator::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef SceneGenerator::Point3DContainer Point3DContainer;
  typedef SceneGenerator::ImageContainer ImageContainer;
  typedef KeysetGenerator::KeysetContainer KeysetContainer;

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

  SceneGenerator scene_generator(
    focal_length_in_metre,
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
    north_west_angle);

  IntrinsicParamsContainer intrinsic_params_set;
  ExtrinsicParamsContainer extrinsic_params_set;
  ImageContainer images;
  Point3DContainer points;
  ASSERT_EQ(0, scene_generator(intrinsic_params_set,
                               extrinsic_params_set,
                               images,
                               points));

  KeysetGenerator keyset_generator(image_width, image_height);
  KeysetContainer keysets;
  hs::sfm::TrackContainer tracks;
  hs::sfm::CameraViewContainer camera_views;
  ASSERT_EQ(0, keyset_generator(intrinsic_params_set,
                                extrinsic_params_set,
                                points,
                                keysets,
                                tracks,
                                camera_views));

  TestMatchesTracksConvertor tester;
  ASSERT_EQ(0, tester.Test(tracks));
}

}
