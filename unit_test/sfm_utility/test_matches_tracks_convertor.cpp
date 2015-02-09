#include <gtest/gtest.h>

#include <algorithm>

#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/synthetic/keyset_generator.hpp"
#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"

namespace
{

class TestMatchesTracksConvertor
{
public:
  typedef int Err;

  Err Test(const hs::sfm::TrackContainer& tracks_true) const
  {
    typedef std::pair<size_t, size_t> View;
    typedef std::map<View, size_t> ViewTrackIndexer;
    //Generate error tracks and expected tracks
    hs::sfm::TrackContainer tracks_error = tracks_true;
    hs::sfm::TrackContainer tracks_expected = tracks_true;

    ViewTrackIndexer view_track_indexer_true;
    auto itr_track_true = tracks_true.begin();
    auto itr_track_true_end = tracks_true.end();
    for (size_t i = 0; itr_track_true != itr_track_true_end;
         ++itr_track_true, ++i)
    {
      auto itr_view_true = itr_track_true->begin();
      auto itr_view_true_end = itr_track_true->end();
      for (; itr_view_true != itr_view_true_end; ++itr_view_true)
      {
        view_track_indexer_true[*itr_view_true] = i;
      }
    }

    const size_t number_of_error_views = 10;
    size_t error_count = 0;
    for (size_t i = 0; i < tracks_true.size(); i++)
    {
      if (!tracks_true[i].empty())
      {
        View view_true = tracks_true[i][0];

        View view_error = view_true;
        view_error.second++;
        auto itr_view_error = view_track_indexer_true.find(view_error);
        if (itr_view_error != view_track_indexer_true.end() &&
            error_count < number_of_error_views)
        {
          size_t j = itr_view_error->second;
          tracks_error[i].push_back(view_error);
          tracks_expected[i].clear();
          tracks_expected[j].clear();
          error_count++;
        }
      }
    }

    hs::sfm::TrackContainer tracks_expected_reordered;

    auto itr_track_expected = tracks_expected.begin();
    auto itr_track_expected_end = tracks_expected.end();
    for (; itr_track_expected != itr_track_expected_end; ++itr_track_expected)
    {
      if (itr_track_expected->size() > 1)
      {
        std::sort(itr_track_expected->begin(), itr_track_expected->end());
        tracks_expected_reordered.push_back(*itr_track_expected);
      }
    }
    std::sort(tracks_expected_reordered.begin(),
              tracks_expected_reordered.end());

    hs::sfm::MatchesTracksConvertor convertor;
    hs::sfm::MatchContainer matches;
    if (convertor(tracks_error, matches) != 0) return -1;
    hs::sfm::TrackContainer tracks_estimate;
    if (convertor(matches, tracks_estimate) != 0) return -1;

    auto itr_track_estimate = tracks_estimate.begin();
    auto itr_track_estimate_end = tracks_estimate.end();
    for (; itr_track_estimate != itr_track_estimate_end; ++itr_track_estimate)
    {
      std::sort(itr_track_estimate->begin(), itr_track_estimate->end());
    }
    std::sort(tracks_estimate.begin(), tracks_estimate.end());

    if (tracks_expected_reordered == tracks_estimate)
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

  typedef std::pair<size_t, size_t> View;
  typedef std::map<View, size_t> ViewTrackIndexer;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 25;
  size_t number_of_cameras_in_strips = 40;
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

  hs::sfm::TrackContainer tracks_tight;
  for (size_t i = 0; i < tracks.size(); i++)
  {
    if (tracks[i].size() > 1)
    {
      tracks_tight.push_back(tracks[i]);
    }
  }

  TestMatchesTracksConvertor tester;
  ASSERT_EQ(0, tester.Test(tracks_tight));
}

}
