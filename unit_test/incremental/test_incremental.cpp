#include <gtest/gtest.h>

#include "data_tester.hpp"
#include "synthetic_data_generator.hpp"
#include "real_data_generator.hpp"

#include "hs_test_utility/test_env/data_path.hpp"

#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"
#include "hs_sfm/incremental/incremental.hpp"
#include "hs_sfm/incremental/gcp_similar_transform_estimator.hpp"
#include "hs_sfm/triangulate/multiple_view_maximum_likelihood_estimator.hpp"
#include "hs_sfm/incremental/bundle_adjustment_gcp_constrained_optimizor.hpp"
#define DEBUG_TMP 0
#if DEBUG_TMP
#include "hs_sfm/sfm_utility/debug_tmp.hpp"
#endif

namespace
{

template <typename _Scalar>
class TestIncremental
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  enum ErrorCode
  {
    NO_ERROR = 0,
    IMAGE_INTRINSIC_NOT_MATCH,
    INCREMENTAL_SFM_FAILED,
    REPROJECTIVE_ERROR_TOO_LARGE,
    SIMILAR_TRANSFORM_ESTIMATION_FAILED,
    GCP_CONSTRAINED_OPTIMIZATION_FAILED,
    GCP_ACCURACY_TOO_BAD,
    CHECK_POINT_ACCURACY_TOO_BAD,
    ESTRINSIC_ACCRURACY_TOO_BAD,
    POINT_ACCURACY_TOO_BAD
  };

private:
  typedef hs::sfm::incremental::IncrementalSFM<Scalar> IncrementalSFM;
  typedef hs::sfm::ObjectIndexMap ObjectIndexMap;
  typedef hs::sfm::ViewInfoIndexer ViewInfoIndexer;
  typedef hs::sfm::incremental::DataTester<Scalar> Tester;
  typedef hs::sfm::incremental::GCPSimilarTransformEstimator<Scalar>
          SimilarTransformEstimator;
  typedef typename SimilarTransformEstimator::Rotation Rotation;
  typedef typename SimilarTransformEstimator::Translate Translate;
  typedef hs::sfm::triangulate::MultipleViewMaximumLikelihoodEstimator<Scalar>
          TriangulateMLEstimator;
  typedef typename TriangulateMLEstimator::Key TriangulateKey;
  typedef typename TriangulateMLEstimator::KeyContainer TriangulateKeyContainer;
  typedef hs::sfm::incremental::BundleAdjustmentGCPConstrainedOptimizor<Scalar>
          GCPConstrainedOptimizor;
#if DEBUG_TMP
  typedef typename hs::sfm::DebugTrue<Scalar> DebugTrueType;
#endif

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
    Scalar ground_resolution,
    Scalar key_stddev,
    const std::string& test_name)
    : key_stddev_(key_stddev),
      test_name_(test_name),
      ground_resolution_(ground_resolution) {}

  Err operator() (
    const IntrinsicParamsContainer& intrinsic_params_set_initial,
    const std::vector<size_t>& image_intrinsic_map_input,
    const MatchContainer& matches,
    const KeysetContainer& keysets,
    const PointContainer& gcps,
    const TrackContainer& tracks_gcp,
    const KeysetContainer& keysets_gcp,
    const PointContainer& check_points,
    const hs::sfm::TrackContainer& tracks_check_point,
    const KeysetContainer& keysets_checkout_point_noised,
    const IntrinsicParamsContainer& intrinsic_params_set_true =
      IntrinsicParamsContainer(),
    const ExtrinsicParamsContainer& extrinsic_params_set_absolute_true =
      ExtrinsicParamsContainer(),
    const PointContainer& points_absolute_true =
      PointContainer()) const
  {
    if (keysets.size() != image_intrinsic_map_input.size())
      return IMAGE_INTRINSIC_NOT_MATCH;
    ExtrinsicParamsContainer extrinsic_params_set_relative_estimate;
    ObjectIndexMap image_intrinsic_map(image_intrinsic_map_input.size());
    for (size_t i = 0; i < image_intrinsic_map_input.size(); i++)
    {
      image_intrinsic_map.SetObjectId(i, image_intrinsic_map_input[i]);
    }
    ObjectIndexMap image_extrinsic_map;
    PointContainer points_relative_estimate;
    TrackContainer tracks;
    ObjectIndexMap track_point_map;
    ViewInfoIndexer view_info_indexer;
    IncrementalSFM incremental_sfm(100, 8, 2, 7);
    IntrinsicParamsContainer intrinsic_params_set_relative_estimate =
      intrinsic_params_set_initial;
#if DEBUG_TMP
    DebugTrueType debug_true(intrinsic_params_set_true,
                             extrinsic_params_set_absolute_true,
                             points_absolute_true,
                             test_name_);
#endif
    if (incremental_sfm(image_intrinsic_map,
                        matches,
                        keysets,
                        intrinsic_params_set_relative_estimate,
                        extrinsic_params_set_relative_estimate,
                        image_extrinsic_map,
                        points_relative_estimate,
                        tracks,
                        track_point_map,
                        view_info_indexer
#if DEBUG_TMP
                        , NULL, debug_true
#endif
                        ) != 0)
    {
      return INCREMENTAL_SFM_FAILED;
    }

    Err result = 0;
    while (1)
    {
      result = TestReprojectiveError(keysets,
                                     intrinsic_params_set_relative_estimate,
                                     extrinsic_params_set_relative_estimate,
                                     points_relative_estimate,
                                     tracks,
                                     image_intrinsic_map,
                                     image_extrinsic_map,
                                     track_point_map,
                                     view_info_indexer);
      if (result != 0) break;

      IntrinsicParamsContainer intrinsic_params_set_absolute_estimate;
      ExtrinsicParamsContainer extrinsic_params_set_absolute_estimate;
      PointContainer points_absolute_estimate;
      result = TestGCPError(keysets,
                            intrinsic_params_set_relative_estimate,
                            extrinsic_params_set_relative_estimate,
                            points_relative_estimate,
                            tracks,
                            image_intrinsic_map,
                            image_extrinsic_map,
                            track_point_map,
                            view_info_indexer,
                            keysets_gcp,
                            tracks_gcp,
                            gcps,
                            intrinsic_params_set_absolute_estimate,
                            extrinsic_params_set_absolute_estimate,
                            points_absolute_estimate);
      if (result != 0) break;

      result = TestCheckPointsError(intrinsic_params_set_absolute_estimate,
                                    extrinsic_params_set_absolute_estimate,
                                    image_intrinsic_map,
                                    image_extrinsic_map,
                                    check_points,
                                    tracks_check_point,
                                    keysets_checkout_point_noised);
      if (result != 0) break;

      result = TestExtrinsicError(extrinsic_params_set_absolute_estimate,
                                  image_extrinsic_map,
                                  extrinsic_params_set_absolute_true);
      if (result != 0) break;

      result = TestPointsError(points_absolute_estimate,
                               track_point_map,
                               points_absolute_true);
      if (result != 0) break;

      break;
    }

    return result;
  }

private:
  Err TestReprojectiveError(
    const KeysetContainer& keysets,
    const IntrinsicParamsContainer& intrinsic_params_set_relative_estimate,
    const ExtrinsicParamsContainer& extrinsic_params_set_relative_estimate,
    const PointContainer& points_relative_estimate,
    const hs::sfm::TrackContainer& tracks,
    const hs::sfm::ObjectIndexMap& image_intrinsic_map,
    const hs::sfm::ObjectIndexMap& image_extrinsic_map,
    const hs::sfm::ObjectIndexMap& track_point_map,
    const hs::sfm::ViewInfoIndexer& view_info_indexer) const
  {
    Tester tester;
    if (tester.TestReprojectiveError(
          keysets,
          intrinsic_params_set_relative_estimate,
          tracks,
          image_intrinsic_map,
          image_extrinsic_map,
          track_point_map,
          view_info_indexer,
          extrinsic_params_set_relative_estimate,
          points_relative_estimate,
          key_stddev_) != 0)
    {
      return REPROJECTIVE_ERROR_TOO_LARGE;
    }
    else
    {
      return 0;
    }
  }

  Err TestGCPError(
    const KeysetContainer& keysets,
    const IntrinsicParamsContainer& intrinsic_params_set_relative_estimate,
    const ExtrinsicParamsContainer& extrinsic_params_set_relative_estimate,
    const PointContainer& points_relative_estimate,
    const hs::sfm::TrackContainer& tracks,
    const hs::sfm::ObjectIndexMap& image_intrinsic_map,
    const hs::sfm::ObjectIndexMap& image_extrinsic_map,
    const hs::sfm::ObjectIndexMap& track_point_map,
    const hs::sfm::ViewInfoIndexer& view_info_indexer,
    const KeysetContainer& keysets_gcp,
    const hs::sfm::TrackContainer& tracks_gcp,
    const PointContainer& gcps,
    IntrinsicParamsContainer& intrinsic_params_set_absolute_estimate,
    ExtrinsicParamsContainer& extrinsic_params_set_absolute_estimate,
    PointContainer& points_absolute_estimate) const
  {
    SimilarTransformEstimator similar_transform_estimator;
    ObjectIndexMap track_point_map_gcp;
    PointContainer gcps_relative;
    Rotation rotation_similar;
    Translate translate_similar;
    Scalar scale_similar;
    if (similar_transform_estimator(keysets_gcp,
                                    intrinsic_params_set_relative_estimate,
                                    extrinsic_params_set_relative_estimate,
                                    tracks_gcp,
                                    gcps,
                                    image_intrinsic_map,
                                    image_extrinsic_map,
                                    4,
                                    key_stddev_ * Scalar(4),
                                    rotation_similar,
                                    translate_similar,
                                    scale_similar,
                                    track_point_map_gcp,
                                    gcps_relative) != 0)
    {
      return SIMILAR_TRANSFORM_ESTIMATION_FAILED;
    }

    intrinsic_params_set_absolute_estimate =
      intrinsic_params_set_relative_estimate;

    extrinsic_params_set_absolute_estimate =
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

    points_absolute_estimate = points_relative_estimate;
    size_t number_of_points = points_absolute_estimate.size();
    for (size_t i = 0; i < number_of_points; i++)
    {
      Point& point = points_absolute_estimate[i];
      point = scale_similar * (rotation_similar * point) + translate_similar;
    }

    PointContainer gcps_absolute_estimate;
    ObjectIndexMap estimate_measure_map;
    GCPConstrainedOptimizor gcp_constrained_optimizor(
                              1, Scalar(0.005), Scalar(0.01));
    if (gcp_constrained_optimizor(keysets,
                                  image_intrinsic_map,
                                  tracks,
                                  image_extrinsic_map,
                                  track_point_map,
                                  view_info_indexer,
                                  keysets_gcp,
                                  tracks_gcp,
                                  gcps,
                                  intrinsic_params_set_absolute_estimate,
                                  extrinsic_params_set_absolute_estimate,
                                  points_absolute_estimate,
                                  gcps_absolute_estimate,
                                  estimate_measure_map) != 0)
    {
      return GCP_CONSTRAINED_OPTIMIZATION_FAILED;
    }

    //size_t number_of_available_gcps = gcps_relative.size();
    //for (size_t i = 0; i < number_of_available_gcps; i++)
    //{
    //  Point& gcp = gcps_absolute_estimate[i];
    //  gcp = scale_similar * (rotation_similar * gcp) + translate_similar;
    //}

    size_t number_of_gcps_estimate = gcps_absolute_estimate.size();
    PointContainer gcps_absolute_true_reordered(number_of_gcps_estimate);
    for (size_t i = 0; i < number_of_gcps_estimate; i++)
    {
      if (estimate_measure_map.IsValid(i))
      {
        gcps_absolute_true_reordered[i] = gcps[estimate_measure_map[i]];
      }
    }

    std::string gcp_accuracy_incremental_path =
      test_name_ + "_gcp_accuracy_incremental.txt";
    Tester tester;
    if (tester.TestPointsAccuracy(
        gcps_absolute_true_reordered,
        gcps_absolute_estimate,
        gcp_accuracy_incremental_path,
        ground_resolution_ * 4) != 0) return GCP_ACCURACY_TOO_BAD;

    return 0;
  }

  Err TestCheckPointsError(
    const IntrinsicParamsContainer& intrinsic_params_set_absolute_estimate,
    const ExtrinsicParamsContainer& extrinsic_params_set_absolute_estimate,
    const hs::sfm::ObjectIndexMap& image_intrinsic_map,
    const hs::sfm::ObjectIndexMap& image_extrinsic_map,
    const PointContainer& check_points,
    const hs::sfm::TrackContainer& tracks_check_point,
    const KeysetContainer& keysets_check_point) const
  {
    PointContainer check_points_true;
    PointContainer check_points_estimate;
    TriangulateMLEstimator triangulate_estimator;
    for (size_t i = 0; i < tracks_check_point.size(); i++)
    {
      size_t track_size = tracks_check_point[i].size();
      TriangulateKeyContainer keys;
      IntrinsicParamsContainer intrinsic_params_set_triangulate;
      ExtrinsicParamsContainer extrinsic_params_set_triangulate;
      for (size_t j = 0; j < track_size; j++)
      {
        size_t image_id = tracks_check_point[i][j].first;
        size_t key_id = tracks_check_point[i][j].second;
        if (!image_intrinsic_map.IsValid(image_id) ||
            !image_extrinsic_map.IsValid(image_id)) continue;
        size_t intrinsic_id = image_intrinsic_map[image_id];
        size_t extrinsic_id = image_extrinsic_map[image_id];
        intrinsic_params_set_triangulate.push_back(
          intrinsic_params_set_absolute_estimate[intrinsic_id]);
        extrinsic_params_set_triangulate.push_back(
          extrinsic_params_set_absolute_estimate[extrinsic_id]);
        keys.push_back(keysets_check_point[image_id][key_id]);
      }
      if (keys.size() > 1)
      {
        Point point;
        if (triangulate_estimator(intrinsic_params_set_triangulate,
                                  extrinsic_params_set_triangulate,
                                  keys,
                                  point) == 0)
        {
          check_points_estimate.push_back(point);
          check_points_true.push_back(check_points[i]);
        }
      }
    }

    std::string check_point_accuracy_incremental_path =
      test_name_ + "_check_point_accuracy_incremental.txt";
    Tester tester;
    if (tester.TestPointsAccuracy(
        check_points_true,
        check_points_estimate,
        check_point_accuracy_incremental_path,
        ground_resolution_ * 4) != 0) return CHECK_POINT_ACCURACY_TOO_BAD;

    return 0;
  }

  Err TestExtrinsicError(
    const ExtrinsicParamsContainer& extrinsic_params_set_absolute_estimate,
    const hs::sfm::ObjectIndexMap& image_extrinsic_map,
    const ExtrinsicParamsContainer& extrinsic_params_set_absolute_true) const
  {
    if (!extrinsic_params_set_absolute_true.empty())
    {
      std::string extrinsic_accuracy_incremental_path =
        test_name_ + "_extrinsic_accuracy_incremental.txt";
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

      Tester tester;
      if (tester.TestExtrinsicAccuracy(
            extrinsic_params_set_absolute_true_reordered,
            extrinsic_params_set_absolute_estimate,
            extrinsic_accuracy_incremental_path,
            ground_resolution_ * 8) != 0)
      {
        return ESTRINSIC_ACCRURACY_TOO_BAD;
      }
      else
      {
        return 0;
      }
    }
    else
    {
      return 0;
    }
  }

  Err TestPointsError(
    const PointContainer& points_absolute_estimate,
    const hs::sfm::ObjectIndexMap& track_point_map,
    const PointContainer& points_absolute_true) const
  {
    if (!points_absolute_true.empty())
    {
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

      std::string point_accuracy_incremental_path =
        test_name_ + "_point_accuracy_incremental.txt";
      Tester tester;
      if (tester.TestPointsAccuracy(
            points_absolute_true_reordered,
            points_absolute_estimate,
            point_accuracy_incremental_path,
            ground_resolution_ * 4) != 0)
      {
        return POINT_ACCURACY_TOO_BAD;
      }
      else
      {
        return 0;
      }
    }
    else
    {
      return 0;
    }
  }

private:
  Scalar key_stddev_;
  std::string test_name_;
  Scalar ground_resolution_;
};

TEST(TestIncremental, Synthetic120ImagesTest)
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

  Tester tester(0.1, key_stddev, test_name);
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

TEST(TestIncremental, Synthetic240ImagesTest)
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

  typedef Generator::FlightGenerator FlightGenerator;
  typedef Generator::FlightGeneratorContainer FlightGeneratorContainer;
  typedef Generator::Image Image;
  typedef Generator::ImageContainer ImageContainer;

  FlightGeneratorContainer flight_generators;
  IntrinsicParamsContainer intrinsic_params_set_true;
  IntrinsicParamsContainer intrinsic_params_set_initial;

  Scalar focal_length_in_metre_0 = 0.018858358970276164;
  size_t number_of_strips_0 = 2;
  size_t number_of_cameras_in_strip_0 = 60;
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
  size_t number_of_cameras_in_strip_1 = 60;
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
  size_t number_of_cameras_in_strip_2 = 60;
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
  size_t number_of_points = 100000;
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

  std::string test_name = "synthetic_data_240_images";

  Tester tester(0.1, key_stddev, test_name);
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

TEST(TestIncremental, Synthetic480ImagesTest)
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

  typedef Generator::FlightGenerator FlightGenerator;
  typedef Generator::FlightGeneratorContainer FlightGeneratorContainer;
  typedef Generator::Image Image;
  typedef Generator::ImageContainer ImageContainer;

  FlightGeneratorContainer flight_generators;
  IntrinsicParamsContainer intrinsic_params_set_true;
  IntrinsicParamsContainer intrinsic_params_set_initial;

  Scalar focal_length_in_metre_0 = 0.018858358970276164;
  size_t number_of_strips_0 = 4;
  size_t number_of_cameras_in_strip_0 = 60;
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
  size_t number_of_strips_1 = 4;
  size_t number_of_cameras_in_strip_1 = 60;
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
  size_t number_of_strips_2 = 4;
  size_t number_of_cameras_in_strip_2 = 60;
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
  size_t number_of_points = 200000;
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

  std::string test_name = "synthetic_data_480_images";

  Tester tester(0.1, key_stddev, test_name);
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

TEST(TestIncremental, Synthetic960ImagesTest)
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

  typedef Generator::FlightGenerator FlightGenerator;
  typedef Generator::FlightGeneratorContainer FlightGeneratorContainer;
  typedef Generator::Image Image;
  typedef Generator::ImageContainer ImageContainer;

  FlightGeneratorContainer flight_generators;
  IntrinsicParamsContainer intrinsic_params_set_true;
  IntrinsicParamsContainer intrinsic_params_set_initial;

  Scalar focal_length_in_metre_0 = 0.018858358970276164;
  size_t number_of_strips_0 = 8;
  size_t number_of_cameras_in_strip_0 = 60;
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
  size_t number_of_strips_1 = 8;
  size_t number_of_cameras_in_strip_1 = 60;
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
  size_t number_of_strips_2 = 8;
  size_t number_of_cameras_in_strip_2 = 60;
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
  size_t number_of_points = 400000;
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

  std::string test_name = "synthetic_data_960_images";

  Tester tester(0.1, key_stddev, test_name);
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

//TEST(TestIncremental, Real120ImageTest)
//{
//  typedef double Scalar;
//  typedef size_t ImageDimension;
//  typedef TestIncremental<Scalar> Tester;
//  typedef Tester::IntrinsicParams IntrinsicParams;
//  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
//  typedef Tester::KeysetContainer KeysetContainer;
//  typedef Tester::PointContainer PointContainer;

//  typedef hs::sfm::incremental::RealDataGenerator<Scalar> Generator;

//  typedef hs::sfm::MatchContainer MatchContainer;
//  typedef hs::sfm::TrackContainer TrackContainer;

//  std::string data_path = hs::test::getTestDataPath();
//  std::string out_path =
//    data_path + "sfm/incremental/real_data_120_images/bundler.out";
//  std::string gcp_path =
//    data_path + "sfm/incremental/real_data_120_images/gcp.xml";
//  IntrinsicParamsContainer intrinsic_params_set_initial;
//  intrinsic_params_set_initial.push_back(IntrinsicParams(7692.30769230769231,
//                                                         0,
//                                                         3000,
//                                                         2000));
//  intrinsic_params_set_initial.push_back(IntrinsicParams(7692.30769230769231,
//                                                         0,
//                                                         3000,
//                                                         2000));
//  intrinsic_params_set_initial.push_back(IntrinsicParams(7692.30769230769231,
//                                                         0,
//                                                         3000,
//                                                         2000));
//  std::vector<size_t> image_intrinsic_map(123);
//  for (size_t i = 0; i < 2; i++)
//  {
//    image_intrinsic_map[i] = 0;
//  }
//  for (size_t i = 2; i < 20; i++)
//  {
//    image_intrinsic_map[i] = 1;
//  }
//  for (size_t i = 20; i < 123; i++)
//  {
//    image_intrinsic_map[i] = 2;
//  }

//  KeysetContainer keysets;
//  MatchContainer matches;
//  PointContainer gcps;
//  TrackContainer tracks_gcp;
//  KeysetContainer keysets_gcp;
//  ASSERT_EQ(0, Generator::LoadBundlerOutFile(out_path, 6000, 4000,
//                                             keysets, matches));
//  ASSERT_EQ(0, Generator::LoadGCPs(gcp_path, gcps, tracks_gcp, keysets_gcp));

//  Scalar key_stddev = Scalar(1);
//  std::string test_name = "real_data";
//  Tester tester(0.1, key_stddev, test_name);
//  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
//                      image_intrinsic_map,
//                      matches,
//                      keysets,
//                      gcps,
//                      tracks_gcp,
//                      keysets_gcp));
//}

}
