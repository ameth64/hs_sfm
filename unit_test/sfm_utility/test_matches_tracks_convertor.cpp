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

  Err Test(const hs::sfm::TrackContainer& tracks_true,
           const hs::sfm::TrackContainer& tracks_error) const
  {
    hs::sfm::MatchesTracksConvertor convertor;
    hs::sfm::MatchContainer matches;
    if (convertor(tracks_error, matches) != 0) return -1;
    hs::sfm::TrackContainer tracks_estimate;
    if (convertor(matches, tracks_estimate) != 0) return -1;

    hs::sfm::TrackContainer tracks_reordered;
    auto itr_track_true = tracks_true.begin();
    auto itr_track_true_end = tracks_true.end();
    for (; itr_track_true != itr_track_true_end; ++itr_track_true)
    {
      //if (itr_track_true->size() > 1)
      {
        hs::sfm::Track track = *itr_track_true;
        std::sort(track.begin(), track.end());
        tracks_reordered.push_back(track);
      }
    }
    std::sort(tracks_reordered.begin(), tracks_reordered.end());

    auto itr_track_estimate = tracks_estimate.begin();
    auto itr_track_estimate_end = tracks_estimate.end();
    for (; itr_track_estimate != itr_track_estimate_end; ++itr_track_estimate)
    {
      std::sort(itr_track_estimate->begin(), itr_track_estimate->end());
    }
    std::sort(tracks_estimate.begin(), tracks_estimate.end());

    for (size_t i = 0; i < std::max(tracks_true.size(),
                                    tracks_error.size()); i++)
    {
      if (i < tracks_true.size())
      {
        std::cout<<"track_true["<<i<<"]:\n";
        for (size_t j = 0; j < tracks_true[i].size(); j++)
        {
          std::cout<<tracks_true[i][j].first<<" "
                   <<tracks_true[i][j].second<<"\n";
        }
      }
      if (i < tracks_error.size())
      {
        std::cout<<"track_error["<<i<<"]:\n";
        for (size_t j = 0; j < tracks_error[i].size(); j++)
        {
          std::cout<<tracks_error[i][j].first<<" "
                   <<tracks_error[i][j].second<<"\n";
        }
      }
    }

    hs::sfm::TrackContainer tracks_error_reordered = tracks_error;
    auto itr_track_error_reordered = tracks_error_reordered.begin();
    auto itr_track_error_reordered_end = tracks_error_reordered.end();
    for (; itr_track_error_reordered != itr_track_error_reordered_end;
         ++itr_track_error_reordered)
    {
      std::sort(itr_track_error_reordered->begin(),
                itr_track_error_reordered->end());
    }
    std::sort(tracks_error_reordered.begin(), tracks_error_reordered.end());

    std::cout<<"tracks_error_reordered.size():"
             <<tracks_error_reordered.size()<<"\n";
    std::cout<<"tracks_estimate.size():"<<tracks_estimate.size()<<"\n";

    for (size_t i = 0; i < std::max(tracks_error_reordered.size(),
                                    tracks_estimate.size()); i++)
    {
      if (i < tracks_error_reordered.size())
      {
        std::cout<<"track_error_reordered["<<i<<"]:\n";
        for (size_t j = 0; j < tracks_error_reordered[i].size(); j++)
        {
          std::cout<<tracks_error_reordered[i][j].first<<" "
                   <<tracks_error_reordered[i][j].second<<"\n";
        }
      }
      if (i < tracks_estimate.size())
      {
        std::cout<<"track_estimate["<<i<<"]:\n";
        for (size_t j = 0; j < tracks_estimate[i].size(); j++)
        {
          std::cout<<tracks_estimate[i][j].first<<" "
                   <<tracks_estimate[i][j].second<<"\n";
        }
      }
    }

    if (tracks_reordered == tracks_estimate)
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
  size_t number_of_strips = 7;
  size_t number_of_cameras_in_strips = 15;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 20;
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

  hs::sfm::TrackContainer tracks_error = tracks_tight;
  hs::sfm::TrackContainer tracks_true = tracks_tight;
  size_t number_of_tracks = tracks_tight.size();
  for (size_t i = 0; i < number_of_tracks; i++)
  {
    size_t j = i + 1;
    if (i % 1000 == 4 && !tracks_tight[i].empty())
    {
      tracks_error[i].push_back(tracks_tight[j][0]);
      tracks_error[j].push_back(tracks_tight[i][0]);
      tracks_true[i].clear();
      tracks_true[j].clear();
    }
    if (i % 1000 == 6 && !tracks_tight[i].empty())
    {
      tracks_error[i].push_back(tracks_tight[i][0]);
      tracks_error[i].rbegin()->second++;
      tracks_true[i].clear();
    }
  }

  TestMatchesTracksConvertor tester;
  ASSERT_EQ(0, tester.Test(tracks_true, tracks_error));
}

}
