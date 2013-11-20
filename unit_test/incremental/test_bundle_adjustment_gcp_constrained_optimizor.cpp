#include <sstream>

#include <gtest/gtest.h>

#include "hs_test_utility/test_env/data_path.hpp"

#include "hs_sfm/sfm_file_io/extrinsic_params_set_loader.hpp"
#include "hs_sfm/sfm_file_io/intrinsic_params_set_loader.hpp"
#include "hs_sfm/sfm_file_io/keyset_loader.hpp"
#include "hs_sfm/sfm_file_io/object_index_loader.hpp"
#include "hs_sfm/sfm_file_io/points_loader.hpp"
#include "hs_sfm/sfm_file_io/tracks_loader.hpp"
#include "hs_sfm/incremental/gcp_similar_transform_estimator.hpp"

#include "data_tester.hpp"

namespace
{

template <typename _Scalar>
class TestBundleAdjustmentGCPConstrainedOptimizor
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

private:
  typedef hs::sfm::incremental::GCPSimilarTransformEstimator<Scalar>
          SimilarTransformEstimator;

  typedef hs::sfm::incremental::DataTester<Scalar> Tester;

public:
  typedef typename SimilarTransformEstimator::ExtrinsicParams ExtrinsicParams;
  typedef typename SimilarTransformEstimator::IntrinsicParams IntrinsicParams;
  typedef typename SimilarTransformEstimator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename SimilarTransformEstimator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename SimilarTransformEstimator::Keyset Keyset;
  typedef typename SimilarTransformEstimator::KeysetContainer KeysetContainer;
  typedef typename SimilarTransformEstimator::Point Point;
  typedef typename SimilarTransformEstimator::PointContainer PointContainer;

private:
  typedef typename SimilarTransformEstimator::Rotation Rotation;
  typedef typename SimilarTransformEstimator::Translate Translate;

public:
  Err operator() (
    const KeysetContainer& keysets_noised,
    const IntrinsicParamsContainer& intrinsic_params_set,
    const ExtrinsicParamsContainer& extrinsic_params_set_absolute_true,
    const PointContainer& points_absolute_true,
    const hs::sfm::TrackContainer& tracks,
    const hs::sfm::ViewInfoIndexer& view_info_indexer,
    const hs::sfm::ObjectIndexMap& image_extrinsic_map,
    const hs::sfm::ObjectIndexMap& track_point_map,
    const ExtrinsicParamsContainer& extrinsic_params_set_relative_estimate,
    const PointContainer& points_relative_estimate,
    const KeysetContainer& keysets_gcp,
    const hs::sfm::TrackContainer& tracks_gcp,
    const PointContainer& gcps,
    Scalar key_stddev,
    const std::string& test_name) const
  {
    SimilarTransformEstimator similar_transform_estimator;
    hs::sfm::ObjectIndexMap track_point_map_gcp;
    hs::sfm::ViewInfoIndexer view_info_indexer_gcp;
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
                                    key_stddev * Scalar(4),
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
      test_name + "_extrinsic_accuracy_similar_transform.txt";

    Tester tester;

    if (tester.TestExtrinsicAccuracy(
          extrinsic_params_set_absolute_true_reordered,
          extrinsic_params_set_absolute_estimate,
          extrinsic_accuracy_similar_transform_path,
          Scalar(2)) != 0) return -1;

    std::string point_accuracy_similar_transform_path =
      test_name + "_point_accuracy_similar_transform.txt";
    if (tester.TestPointsAccuracy(
          points_absolute_true_reordered,
          points_absolute_estimate,
          point_accuracy_similar_transform_path,
          Scalar(3)) != 0) return -1;

    return 0;
  }
};

TEST(TestBundleAdjustmentGCPConstrainedOptimizor, BigDataTest)
{
  typedef double Scalar;
  typedef int Err;

  typedef TestBundleAdjustmentGCPConstrainedOptimizor<Scalar> Tester;
  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Tester::KeysetContainer KeysetContainer;
  typedef Tester::PointContainer PointContainer;

  std::string test_name = "big_data";
  std::string test_prefix = hs::test::getTestDataPath();
  test_prefix += "sfm/incremental/" + test_name;

  hs::sfm::fileio::ExtrinsicParamsSetLoader<Scalar> extrinsic_set_loader;
  hs::sfm::fileio::PointsLoader<Scalar> points_loader;
  hs::sfm::fileio::TracksLoader tracks_loader;
  hs::sfm::fileio::IntrinsicParamsSetLoader<Scalar> intrinsic_set_loader;
  hs::sfm::fileio::ObjectIndexLoader object_index_loader;
  hs::sfm::fileio::KeysetLoader<Scalar> keyset_loader;

  std::string extrinsic_set_absolute_true_path =
    test_prefix + "_extrinsic_params_set_absolute_true.txt";
  ExtrinsicParamsContainer extrinsic_params_set_absolute_true;
  ASSERT_EQ(0, extrinsic_set_loader(extrinsic_set_absolute_true_path,
                                    extrinsic_params_set_absolute_true));

  std::string points_absolute_true_path =
    test_prefix + "_points_absolute_true.txt";
  PointContainer points_absolute_true;
  ASSERT_EQ(0, points_loader(points_absolute_true_path,
                             points_absolute_true));

  std::string gcps_path = test_prefix + "_gcps.txt";
  PointContainer gcps;
  ASSERT_EQ(0, points_loader(gcps_path,
                             gcps));

  size_t number_of_images = extrinsic_params_set_absolute_true.size();
  KeysetContainer keysets_noised(number_of_images);
  for (size_t i = 0; i < number_of_images; i++)
  {
    std::stringstream ss;
    ss<<test_prefix<<"_keyset_"<<i<<".txt";
    std::string keyset_path;
    ss>>keyset_path;
    keyset_loader(keyset_path, keysets_noised[i]);
  }

  KeysetContainer keysets_gcp_noised(number_of_images);
  for (size_t i = 0; i < number_of_images; i++)
  {
    std::stringstream ss;
    ss<<test_prefix<<"_keyset_gcp_"<<i<<".txt";
    std::string keyset_gcp_path;
    ss>>keyset_gcp_path;
    keyset_loader(keyset_gcp_path, keysets_gcp_noised[i]);
  }

  hs::sfm::ViewInfoIndexer view_info_indexer_gcp;
  hs::sfm::TrackContainer tracks_gcp;
  std::string tracks_gcp_path = test_prefix + "_tracks_gcp.txt";
  ASSERT_EQ(0, tracks_loader(tracks_gcp_path, tracks_gcp,
                             view_info_indexer_gcp));

  hs::sfm::ViewInfoIndexer view_info_indexer;
  hs::sfm::TrackContainer tracks;
  std::string tracks_path = test_prefix + "_tracks.txt";
  ASSERT_EQ(0, tracks_loader(tracks_path, tracks, view_info_indexer));

  IntrinsicParamsContainer intrinsic_params_set;
  std::string intrinsic_set_path =
    test_prefix + "_intrinsic_params_set.txt";
  ASSERT_EQ(0, intrinsic_set_loader(intrinsic_set_path,
                                    intrinsic_params_set));

  std::string extrinsic_set_relative_estimate_path =
    test_prefix + "_extrinsic_params_set_relative_estimate.txt";
  ExtrinsicParamsContainer extrinsic_params_set_relative_estimate;
  ASSERT_EQ(0, extrinsic_set_loader(
                 extrinsic_set_relative_estimate_path,
                 extrinsic_params_set_relative_estimate));

  std::string points_relative_estimate_path =
    test_prefix + "_points_relative_estimate.txt";
  PointContainer points_relative_estimate;
  ASSERT_EQ(0, points_loader(points_relative_estimate_path,
                            points_relative_estimate));

  std::string track_point_map_path =
    test_prefix + "_track_point_map.txt";
  hs::sfm::ObjectIndexMap track_point_map;
  ASSERT_EQ(0, object_index_loader(track_point_map_path,
                                   track_point_map));

  std::string image_extrinsic_map_path =
    test_prefix + "_image_extrinsic_map.txt";
  hs::sfm::ObjectIndexMap image_extrinsic_map;
  ASSERT_EQ(0, object_index_loader(image_extrinsic_map_path,
                                   image_extrinsic_map));

  Tester tester;
  ASSERT_EQ(0, tester(keysets_noised,
                      intrinsic_params_set,
                      extrinsic_params_set_absolute_true,
                      points_absolute_true,
                      tracks,
                      view_info_indexer,
                      image_extrinsic_map,
                      track_point_map,
                      extrinsic_params_set_relative_estimate,
                      points_relative_estimate,
                      keysets_gcp_noised,
                      tracks_gcp,
                      gcps,
                      1,
                      test_name));

}

}
