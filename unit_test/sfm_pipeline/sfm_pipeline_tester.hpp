#ifndef _HS_SFM_UNIT_TEST_SFM_PIPELINE_SFM_PIPELINE_TESTER_HPP_
#define _HS_SFM_UNIT_TEST_SFM_PIPELINE_SFM_PIPELINE_TESTER_HPP_

#include "data_tester.hpp"
#include "synthetic_data_generator.hpp"
#include "real_data_generator.hpp"
#include "real_synthetic_data_generator.hpp"

#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"
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

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _SFMPipeline>
class SFMPipelineTester
{
public:
  typedef _SFMPipeline SFMPipeline;
  typedef typename SFMPipeline::Scalar Scalar;
  typedef int Err;

  enum ErrorCode
  {
    NO_ERROR = 0,
    IMAGE_INTRINSIC_NOT_MATCH,
    SFM_FAILED,
    REPROJECTIVE_ERROR_TOO_LARGE,
    SIMILAR_TRANSFORM_ESTIMATION_FAILED,
    GCP_CONSTRAINED_OPTIMIZATION_FAILED,
    GCP_ACCURACY_TOO_BAD,
    CHECK_POINT_ACCURACY_TOO_BAD,
    EXTRINSIC_ACCRURACY_TOO_BAD,
    POINT_ACCURACY_TOO_BAD
  };

private:
  typedef hs::sfm::ObjectIndexMap ObjectIndexMap;
  typedef hs::sfm::ViewInfoIndexer ViewInfoIndexer;
  typedef hs::sfm::pipeline::DataTester<Scalar> Tester;
  typedef hs::sfm::pipeline::GCPSimilarTransformEstimator<Scalar>
          SimilarTransformEstimator;
  typedef typename SimilarTransformEstimator::Rotation Rotation;
  typedef typename SimilarTransformEstimator::Translate Translate;
  typedef hs::sfm::triangulate::MultipleViewMaximumLikelihoodEstimator<Scalar>
          TriangulateMLEstimator;
  typedef typename TriangulateMLEstimator::Key TriangulateKey;
  typedef typename TriangulateMLEstimator::KeyContainer TriangulateKeyContainer;
  typedef hs::sfm::pipeline::BundleAdjustmentGCPConstrainedOptimizor<Scalar>
          GCPConstrainedOptimizor;

public:
  typedef typename SFMPipeline::IntrinsicParams IntrinsicParams;
  typedef typename SFMPipeline::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename SFMPipeline::ExtrinsicParams ExtrinsicParams;
  typedef typename SFMPipeline::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename SFMPipeline::KeysetContainer KeysetContainer;
  typedef typename SFMPipeline::Point Point;
  typedef typename SFMPipeline::PointContainer PointContainer;
  typedef hs::sfm::MatchContainer MatchContainer;
  typedef hs::sfm::TrackContainer TrackContainer;

public:
  SFMPipelineTester(
    SFMPipeline& sfm_pipeline,
    Scalar ground_resolution,
    Scalar key_stddev,
    const std::string& test_name)
    : sfm_pipeline_(sfm_pipeline),
      key_stddev_(key_stddev),
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
    IntrinsicParamsContainer intrinsic_params_set_relative_estimate =
      intrinsic_params_set_initial;
    if (sfm_pipeline_(image_intrinsic_map,
                      matches,
                      keysets,
                      intrinsic_params_set_relative_estimate,
                      extrinsic_params_set_relative_estimate,
                      image_extrinsic_map,
                      points_relative_estimate,
                      tracks,
                      track_point_map,
                      view_info_indexer
                      ) != 0)
    {
      return SFM_FAILED;
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

    PointContainer gcps_absolute_estimate = gcps_relative;
    for (size_t i = 0; i < gcps_absolute_estimate.size(); i++)
    {
      Point& gcp = gcps_absolute_estimate[i];
      gcp = scale_similar * (rotation_similar * gcp) + translate_similar;
    }
    ObjectIndexMap estimate_measure_map = track_point_map_gcp;
    GCPConstrainedOptimizor gcp_constrained_optimizor(
                              1, Scalar(0.005), Scalar(0.005));
    Tester tester;
    if (tester.TestReprojectiveError(
          keysets,
          intrinsic_params_set_absolute_estimate,
          tracks,
          image_intrinsic_map,
          image_extrinsic_map,
          track_point_map,
          view_info_indexer,
          extrinsic_params_set_absolute_estimate,
          points_absolute_estimate,
          key_stddev_) != 0)
    {
      return REPROJECTIVE_ERROR_TOO_LARGE;
    }
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
    if (tester.TestReprojectiveError(
          keysets,
          intrinsic_params_set_absolute_estimate,
          tracks,
          image_intrinsic_map,
          image_extrinsic_map,
          track_point_map,
          view_info_indexer,
          extrinsic_params_set_absolute_estimate,
          points_absolute_estimate,
          key_stddev_) != 0)
    {
      return REPROJECTIVE_ERROR_TOO_LARGE;
    }

#if DEBUG_TMP
    typedef hs::sfm::fileio::ScenePLYSaver<Scalar, size_t> SceneSaver;
    typedef typename SceneSaver::Image Image;
    typedef typename SceneSaver::ImageContainer ImageContainer;
    std::string abs_scene_path = test_name_ + "_abs_scene.ply";
    SceneSaver saver(Scalar(10));
    IntrinsicParamsContainer intrinsic_params_set_out(
      extrinsic_params_set_absolute_estimate.size(),
      intrinsic_params_set_absolute_estimate[0]);
    Image image_out;
    image_out.m_width = 6000;
    image_out.m_height = 4000;
    ImageContainer images_out(extrinsic_params_set_absolute_estimate.size(),
                              image_out);
    saver(abs_scene_path,
          intrinsic_params_set_out,
          extrinsic_params_set_absolute_estimate,
          images_out,
          points_absolute_estimate);

    std::string abs_intrinsic_path = test_name_ + "_abs_intrinsic.txt";
    std::ofstream abs_intrinsic_file(abs_intrinsic_path);
    abs_intrinsic_file.setf(std::ios::fixed);
    abs_intrinsic_file<<std::setprecision(8);
    for (size_t i = 0; i < intrinsic_params_set_absolute_estimate.size(); i++)
    {
      const IntrinsicParams& intrinsic_params =
        intrinsic_params_set_absolute_estimate[i];
      abs_intrinsic_file<<intrinsic_params.focal_length()<<" "
                        <<intrinsic_params.skew()<<" "
                        <<intrinsic_params.principal_point_x()<<" "
                        <<intrinsic_params.principal_point_y()<<" "
                        <<intrinsic_params.pixel_ratio()<<" "
                        <<intrinsic_params.k1()<<" "
                        <<intrinsic_params.k2()<<" "
                        <<intrinsic_params.k3()<<" "
                        <<intrinsic_params.d1()<<" "
                        <<intrinsic_params.d2()<<"\n";
    }
    abs_intrinsic_file.close();

#endif

    size_t number_of_gcps_estimate = gcps_absolute_estimate.size();
    PointContainer gcps_absolute_true_reordered(number_of_gcps_estimate);
    for (size_t i = 0; i < number_of_gcps_estimate; i++)
    {
      if (estimate_measure_map.IsValid(i))
      {
        gcps_absolute_true_reordered[i] = gcps[estimate_measure_map[i]];
      }
    }

    std::string gcp_accuracy_path =
      test_name_ + "_gcp_accuracy.txt";
    if (tester.TestPointsAccuracy(
        gcps_absolute_true_reordered,
        gcps_absolute_estimate,
        gcp_accuracy_path,
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

    std::string check_point_accuracy_path =
      test_name_ + "_check_point_accuracy.txt";
    Tester tester;
    if (tester.TestPointsAccuracy(
        check_points_true,
        check_points_estimate,
        check_point_accuracy_path,
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
      std::string extrinsic_accuracy_path =
        test_name_ + "_extrinsic_accuracy.txt";
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
            extrinsic_accuracy_path,
            ground_resolution_ * 8) != 0)
      {
        return EXTRINSIC_ACCRURACY_TOO_BAD;
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

      std::string point_accuracy_path =
        test_name_ + "_point_accuracy.txt";
      Tester tester;
      if (tester.TestPointsAccuracy(
            points_absolute_true_reordered,
            points_absolute_estimate,
            point_accuracy_path,
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
  SFMPipeline& sfm_pipeline_;
};

}
}
}

#endif
