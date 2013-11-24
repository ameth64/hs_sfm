#include <gtest/gtest.h>

#include "data_tester.hpp"
#include "synthetic_data_generator.hpp"

#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"
#include "hs_sfm/incremental/incremental.hpp"
#include "hs_sfm/incremental/gcp_similar_transform_estimator.hpp"

namespace
{

template <typename _Scalar>
class TestIncremental
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

private:
  typedef hs::sfm::incremental::IncrementalSFM<Scalar> IncrementalSFM;
  typedef hs::sfm::ObjectIndexMap ObjectIndexMap;
  typedef hs::sfm::ViewInfoIndexer ViewInfoIndexer;
  typedef hs::sfm::incremental::DataTester<Scalar> Tester;
  typedef hs::sfm::incremental::GCPSimilarTransformEstimator<Scalar>
          SimilarTransformEstimator;
  typedef typename SimilarTransformEstimator::Rotation Rotation;
  typedef typename SimilarTransformEstimator::Translate Translate;

public:
  typedef typename IncrementalSFM::IntrinsicParams IntrinsicParams;
  typedef typename IncrementalSFM::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename IncrementalSFM::ExtrinsicParams ExtrinsicParams;
  typedef typename IncrementalSFM::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename IncrementalSFM::KeysetContainer KeysetContainer;
  typedef typename IncrementalSFM::Point Point;
  typedef typename IncrementalSFM::PointContainer PointContainer;
  typedef hs::sfm::MatchContainer MatchContainer;
  typedef hs::sfm::TrackContainer TrackContainer;

public:
  TestIncremental(
    Scalar key_stddev,
    const std::string& test_name)
    : key_stddev_(key_stddev),
      test_name_(test_name) {}

  Err operator() (
    const IntrinsicParamsContainer& intrinsic_params_set,
    const MatchContainer& matches,
    const KeysetContainer& keysets,
    const PointContainer& gcps,
    const TrackContainer& tracks_gcp,
    const KeysetContainer& keysets_gcp,
    const ExtrinsicParamsContainer& extrinsic_params_set_absolute_true,
    const PointContainer& points_absolute_true) const
  {

    ExtrinsicParamsContainer extrinsic_params_set_relative_estimate;
    ObjectIndexMap image_extrinsic_map;
    PointContainer points_relative_estimate;
    ObjectIndexMap track_point_map;
    ViewInfoIndexer view_info_indexer;
    IncrementalSFM incremental_sfm;
    if (incremental_sfm(intrinsic_params_set,
                        matches,
                        keysets,extrinsic_params_set_relative_estimate,
                        image_extrinsic_map,
                        points_relative_estimate,
                        track_point_map,
                        view_info_indexer) != 0)
    {
      return -1;
    }

    Tester tester;
    hs::sfm::MatchesTracksConvertor matches_tracks_convertor;
    TrackContainer tracks;
    if (matches_tracks_convertor(matches, tracks) != 0) return -1;
    if (tester.TestReprojectiveError(
          keysets,
          intrinsic_params_set,
          tracks,
          image_extrinsic_map,
          track_point_map,
          view_info_indexer,
          extrinsic_params_set_relative_estimate,
          points_relative_estimate,
          key_stddev_) != 0) return -1;

    SimilarTransformEstimator similar_transform_estimator;
    ObjectIndexMap track_point_map_gcp;
    ViewInfoIndexer view_info_indexer_gcp;
    PointContainer gcps_relative;
    Rotation rotation_similar;
    Translate translate_similar;
    Scalar scale_similar;
    if (similar_transform_estimator(keysets_gcp,
                                    intrinsic_params_set,
                                    extrinsic_params_set_relative_estimate,
                                    tracks_gcp,
                                    gcps,
                                    image_extrinsic_map,
                                    4,
                                    key_stddev_ * Scalar(4),
                                    rotation_similar,
                                    translate_similar,
                                    scale_similar,
                                    track_point_map_gcp,
                                    view_info_indexer_gcp,
                                    gcps_relative) != 0)
    {
      return -1;
    }

    ExtrinsicParamsContainer extrinsic_params_set_absolute_estimate =
      extrinsic_params_set_relative_estimate;
    size_t number_of_extrinsics =
      extrinsic_params_set_absolute_estimate.size();
    for (size_t i = 0; i < number_of_extrinsics; i++)
    {
      ExtrinsicParams& extrinsic_params =
        extrinsic_params_set_absolute_estimate[i];
      extrinsic_params.rotation() =
        extrinsic_params.rotation() * rotation_similar.Inverse();
      extrinsic_params.position() =
        scale_similar * (rotation_similar * extrinsic_params.position()) +
        translate_similar;
    }

    PointContainer points_absolute_estimate = points_relative_estimate;
    size_t number_of_points = points_absolute_estimate.size();
    for (size_t i = 0; i < number_of_points; i++)
    {
      Point& point = points_absolute_estimate[i];
      point = scale_similar * (rotation_similar * point) + translate_similar;
    }

    ExtrinsicParamsContainer extrinsic_params_set_absolute_true_reordered(
                                extrinsic_params_set_absolute_estimate.size());
    size_t number_of_images = extrinsic_params_set_absolute_true.size();
    for (size_t i = 0; i < number_of_images; i++)
    {
      if (image_extrinsic_map.IsValid(i))
      {
        extrinsic_params_set_absolute_true_reordered[image_extrinsic_map[i]] =
          extrinsic_params_set_absolute_true[i];
      }
    }

    PointContainer points_absolute_true_reordered(
                     points_absolute_estimate.size());
    size_t number_of_tracks = points_absolute_true.size();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      if (track_point_map.IsValid(i))
      {
        points_absolute_true_reordered[track_point_map[i]] =
          points_absolute_true[i];
      }
    }

    std::string extrinsic_accuracy_similar_transform_path =
      test_name_ + "_extrinsic_accuracy_incremental.txt";

    if (tester.TestExtrinsicAccuracy(
          extrinsic_params_set_absolute_true_reordered,
          extrinsic_params_set_absolute_estimate,
          extrinsic_accuracy_similar_transform_path,
          Scalar(2)) != 0) return -1;

    std::string point_accuracy_similar_transform_path =
      test_name_ + "_point_accuracy_incremental.txt";
    if (tester.TestPointsAccuracy(
          points_absolute_true_reordered,
          points_absolute_estimate,
          point_accuracy_similar_transform_path,
          Scalar(3)) != 0) return -1;

    return 0;
  }

private:
  Scalar key_stddev_;
  std::string test_name_;
};

TEST(TestIncremental, SyntheticTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestIncremental<Scalar> Tester;
  typedef Tester::IntrinsicParams IntrinsicParams;
  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Tester::ExtrinsicParams ExtrinsicParams;
  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef Tester::KeysetContainer KeysetContainer;
  typedef Tester::PointContainer PointContainer;
  typedef Tester::MatchContainer MatchContainer;
  typedef Tester::TrackContainer TrackContainer;

  typedef hs::sfm::incremental::SyntheticDataGenerator<Scalar, ImageDimension>
          Generator;

  typedef Generator::Image Image;
  typedef Generator::ImageContainer ImageContainer;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 4;
  size_t number_of_cameras_in_strips = 5;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 5000;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 100;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 5;
  Scalar north_west_angle = 60;
  Scalar outlier_ratio = 0.0;
  Scalar key_stddev = 1;

  size_t number_of_gcps = 10;

  Generator generator(focal_length_in_metre,
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
                      key_stddev);

  IntrinsicParamsContainer intrinsic_params_set;
  ExtrinsicParamsContainer extrinsic_params_set_absolute_true;
  ImageContainer images;
  PointContainer points_absolute_true;
  KeysetContainer keysets_true;
  TrackContainer tracks;
  hs::sfm::CameraViewContainer camera_views;

  ASSERT_EQ(0, generator.GenerateAbsoluteScene(
                           intrinsic_params_set,
                           extrinsic_params_set_absolute_true,
                           images,
                           points_absolute_true,
                           keysets_true,
                           tracks,
                           camera_views));

  hs::sfm::MatchesTracksConvertor matches_tracks_convertor;
  MatchContainer matches;
  ASSERT_EQ(0, matches_tracks_convertor(tracks, matches));

  KeysetContainer keysets_noised;
  ASSERT_EQ(0, generator.GenerateNoisedKeysets(keysets_true,
                                               keysets_noised));

  KeysetContainer keysets_gcp_true;
  PointContainer gcps;
  TrackContainer tracks_gcp;
  ASSERT_EQ(0, generator.GenerateGCPData(intrinsic_params_set,
                                         extrinsic_params_set_absolute_true,
                                         number_of_gcps,
                                         keysets_gcp_true,
                                         gcps,
                                         tracks_gcp));

  KeysetContainer keysets_gcp_noised;
  ASSERT_EQ(0, generator.GenerateNoisedKeysets(keysets_gcp_true,
                                               keysets_gcp_noised));

  std::string test_name = "synthetic_data";

  Tester tester(key_stddev, test_name);
  ASSERT_EQ(0, tester(intrinsic_params_set,
                      matches,
                      keysets_noised,
                      gcps,
                      tracks_gcp,
                      keysets_gcp_noised,
                      extrinsic_params_set_absolute_true,
                      points_absolute_true));

}

}
