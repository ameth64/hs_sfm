#ifndef _HS_SFM_SFM_PIPELINE_BUNDLE_ADJUSTMENT_GCP_CONSTRAINED_OPTIMIZOR_HPP_
#define _HS_SFM_SFM_PIPELINE_BUNDLE_ADJUSTMENT_GCP_CONSTRAINED_OPTIMIZOR_HPP_

#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_ceres_optimizor.hpp"
#include "hs_sfm/triangulate/multiple_view_maximum_likelihood_estimator.hpp"

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
class BundleAdjustmentGCPConstrainedOptimizor
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef ImageKeys<Scalar> ImageKeyset;
  typedef EIGEN_STD_VECTOR(ImageKeyset) ImageKeysetContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

  typedef ObjectIndexMap TrackPointMap;
  typedef ObjectIndexMap ImageIntrinsicMap;
  typedef ObjectIndexMap ImageExtrinsicMap;

private:
  typedef hs::sfm::ba::CameraSharedVectorFunction<Scalar> BAVectorFunction;
  typedef typename BAVectorFunction::Index Index;
  typedef typename BAVectorFunction::XVector XVector;
  typedef typename BAVectorFunction::YVector YVector;
  typedef typename BAVectorFunction::FeatureMap FeatureMap;
  typedef typename BAVectorFunction::FeatureMapContainer FeatureMapContainer;
  typedef typename BAVectorFunction::ImageCameraMap ImageCameraMap;
  typedef typename BAVectorFunction::Point BAPoint;
  typedef typename BAVectorFunction::Image BAImage;
  typedef typename BAVectorFunction::Camera BACamera;

  typedef hs::sfm::ba::PointConstraint PointConstraint;
  typedef hs::sfm::ba::PointConstraintContainer PointConstraintContainer;

  typedef hs::sfm::ba::CameraSharedCeresOptimizor<BAVectorFunction>
          BAOptimizor;
  typedef typename BAOptimizor::YCovarianceInverse YCovarianceInverse;
  typedef typename YCovarianceInverse::KeyBlock KeyBlock;

  typedef typename IntrinsicParams::KMatrix KMatrix;
  typedef typename ExtrinsicParams::Position Position;

  typedef hs::sfm::triangulate::MultipleViewMaximumLikelihoodEstimator<Scalar>
          TriangulateMLEstimator;
  typedef typename TriangulateMLEstimator::Key TriangulateKey;
  typedef typename TriangulateMLEstimator::KeyContainer TriangulateKeyContainer;

public:
  BundleAdjustmentGCPConstrainedOptimizor(
    size_t number_of_threads,
    Scalar gcp_planar_accuracy,
    Scalar gcp_height_accuracy,
    Scalar image_accuracy = Scalar(0.1),
    Scalar gcp_image_accuracy = Scalar(4))
    : number_of_threads_(number_of_threads)
    , gcp_planar_accuracy_(gcp_planar_accuracy)
    , gcp_height_accuracy_(gcp_height_accuracy)
    , image_accuracy_(image_accuracy)
    , gcp_image_accuracy_(gcp_image_accuracy) {}

  Err operator() (const ImageKeysetContainer& image_keysets,
                  const ImageIntrinsicMap& image_intrinsic_map,
                  const TrackContainer& tracks,
                  const ImageExtrinsicMap& image_extrinsic_map,
                  const TrackPointMap& track_point_map,
                  const ViewInfoIndexer& view_info_indexer,
                  const ImageKeysetContainer& image_keysets_gcp,
                  const TrackContainer& tracks_gcp,
                  const PointContainer& gcps_measure,
                  IntrinsicParamsContainer& intrinsic_params_set,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  PointContainer& points,
                  PointContainer& gcps_estimate,
                  ObjectIndexMap& estimate_measure_map) const
  {
    size_t number_of_tracks = tracks.size();
    TrackContainer tracks_bundle;
    ObjectIndexMap track_point_map_bundle;
    Index number_of_keys = 0;
    for (size_t track_id = 0; track_id < number_of_tracks; track_id++)
    {
      if (track_point_map.IsValid(track_id))
      {
        size_t number_of_views = tracks[track_id].size();
        hs::sfm::Track track_bundle;
        for (size_t view_id = 0; view_id < number_of_views; view_id++)
        {
          size_t image_id = tracks[track_id][view_id].first;
          const ViewInfo& view_info =
            view_info_indexer.GetViewInfoByTrackImage(track_id, image_id);
          if (image_extrinsic_map.IsValid(image_id) &&
              !view_info.is_blunder)
          {
            track_bundle.push_back(tracks[track_id][view_id]);
          }
        }
        if (!track_bundle.empty())
        {
          tracks_bundle.push_back(track_bundle);
          size_t point_id = track_point_map[track_id];
          track_point_map_bundle.AddObject(point_id);
          number_of_keys += Index(track_bundle.size());
        }
      }
    }

    size_t number_of_gcps = tracks_gcp.size();
    TrackContainer tracks_gcp_estimate;
    TriangulateMLEstimator triangulate_estimator;
    gcps_estimate.clear();
    for (size_t track_id = 0; track_id < number_of_gcps; track_id++)
    {
      size_t number_of_views_gcp = tracks_gcp[track_id].size();

      //三角化像控点作为初始值
      size_t track_size = tracks_gcp[track_id].size();
      TriangulateKeyContainer keys;
      IntrinsicParamsContainer intrinsic_params_set_triangulate;
      ExtrinsicParamsContainer extrinsic_params_set_triangulate;
      hs::sfm::Track track_gcp_estimate;
      for (size_t view_id = 0; view_id < track_size; view_id++)
      {
        size_t image_id = tracks_gcp[track_id][view_id].first;
        size_t key_id = tracks_gcp[track_id][view_id].second;
        if (!image_intrinsic_map.IsValid(image_id) ||
            !image_extrinsic_map.IsValid(image_id)) continue;
        size_t intrinsic_id = image_intrinsic_map[image_id];
        size_t extrinsic_id = image_extrinsic_map[image_id];
        intrinsic_params_set_triangulate.push_back(
          intrinsic_params_set[intrinsic_id]);
        extrinsic_params_set_triangulate.push_back(
          extrinsic_params_set[extrinsic_id]);
        keys.push_back(image_keysets_gcp[image_id][key_id]);
        track_gcp_estimate.push_back(tracks_gcp[track_id][view_id]);
      }
      Point point;
      if (triangulate_estimator(intrinsic_params_set_triangulate,
                                extrinsic_params_set_triangulate,
                                keys,
                                point) != 0) continue;
      gcps_estimate.push_back(point);
      number_of_keys += keys.size();
      estimate_measure_map.AddObject(track_id);
      tracks_gcp_estimate.push_back(track_gcp_estimate);
    }
    size_t number_of_gcps_estimate = gcps_estimate.size();

    //设置feature_maps
    size_t number_of_tracks_bundle = tracks_bundle.size();
    FeatureMapContainer feature_maps;
    for (size_t track_id = 0; track_id < number_of_tracks_bundle; track_id++)
    {
      size_t number_of_views_bundle = tracks_bundle[track_id].size();
      for (size_t view_id = 0; view_id < number_of_views_bundle; view_id++)
      {
        size_t image_id = tracks_bundle[track_id][view_id].first;
        size_t key_id = tracks_bundle[track_id][view_id].second;
        size_t extrinsic_id = image_extrinsic_map[image_id];

        feature_maps.push_back(std::make_pair(extrinsic_id, track_id));
      }
    }
    size_t number_of_feature_map_gcp = 0;
    for (size_t gcp_estimate_id = 0; gcp_estimate_id < gcps_estimate.size();
         gcp_estimate_id++)
    {
      size_t number_of_views = tracks_gcp_estimate[gcp_estimate_id].size();
      for (size_t view_id = 0; view_id < number_of_views; view_id++)
      {
        size_t image_id = tracks_gcp_estimate[gcp_estimate_id][view_id].first;
        size_t key_id = tracks_gcp_estimate[gcp_estimate_id][view_id].second;
        size_t extrinsic_id = image_extrinsic_map[image_id];
        feature_maps.push_back(
          std::make_pair(extrinsic_id,
                         number_of_tracks_bundle + gcp_estimate_id));
        number_of_feature_map_gcp++;
      }
    }
    //设置image_camera_map
    ImageCameraMap image_camera_map;
    ObjectIndexMap extrinsic_image_map(extrinsic_params_set.size());
    for (size_t image_id = 0; image_id < image_extrinsic_map.Size(); image_id++)
    {
      if (image_extrinsic_map.IsValid(image_id))
      {
        extrinsic_image_map[image_extrinsic_map[image_id]] = image_id;
      }
    }
    for (size_t extrinsic_id = 0; extrinsic_id < extrinsic_image_map.Size();
         extrinsic_id++)
    {
      if (!extrinsic_image_map.IsValid(extrinsic_id)) return -1;
      size_t image_id = extrinsic_image_map[extrinsic_id];
      image_camera_map.push_back(image_intrinsic_map[image_id]);
    }

    Index number_of_points = Index(number_of_tracks_bundle +
                                   number_of_gcps_estimate);
    Index number_of_images = Index(extrinsic_params_set.size());
    Index number_of_cameras = Index(intrinsic_params_set.size());
    BAVectorFunction vector_function;
    vector_function.set_number_of_images(number_of_images);
    vector_function.set_number_of_points(number_of_points);
    vector_function.set_number_of_keys(number_of_keys);
    vector_function.set_number_of_cameras(number_of_cameras);
    vector_function.set_feature_maps(feature_maps);
    vector_function.set_image_camera_map(image_camera_map);
    vector_function.intrinsic_computations_mask().set();

    //设置点约束
    PointConstraintContainer point_constraints;
    for (size_t gcp_estimate_id = 0; gcp_estimate_id < number_of_gcps_estimate;
         gcp_estimate_id++)
    {
      PointConstraint point_constraint;
      point_constraint.point_id = number_of_tracks_bundle + gcp_estimate_id;
      point_constraint.mask.set(hs::sfm::ba::POINT_CONSTRAIN_X);
      point_constraint.mask.set(hs::sfm::ba::POINT_CONSTRAIN_Y);
      point_constraint.mask.set(hs::sfm::ba::POINT_CONSTRAIN_Z);
      point_constraints.push_back(point_constraint);
    }
    vector_function.point_constraints() = point_constraints;

    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    XVector initial_x(x_size);
    for (Index track_bundle_id = 0;
         track_bundle_id < Index(number_of_tracks_bundle);
         track_bundle_id++)
    {
      size_t point_id = track_point_map_bundle[track_bundle_id];
      initial_x.template segment<BAVectorFunction::params_per_point_>(
        track_bundle_id * BAVectorFunction::params_per_point_) =
        points[point_id];
    }
    for (size_t gcp_estimate_id = 0; gcp_estimate_id < number_of_gcps_estimate;
         gcp_estimate_id++)
    {
      initial_x.template segment<BAVectorFunction::params_per_point_>(
        Index(number_of_tracks_bundle + gcp_estimate_id) *
        BAVectorFunction::params_per_point_) = gcps_estimate[gcp_estimate_id];
    }

    Index extrinsic_begin = vector_function.GetPointParamsSize();
    for (Index image_id = 0; image_id < number_of_images; image_id++)
    {
      const ExtrinsicParams& extrinsic_params = extrinsic_params_set[image_id];
      Point t = -(extrinsic_params.rotation() * extrinsic_params.position());
      Index offset = extrinsic_begin +
                     image_id * BAVectorFunction::extrinsic_params_per_image_;
      initial_x[offset + 0] = extrinsic_params.rotation()[0];
      initial_x[offset + 1] = extrinsic_params.rotation()[1];
      initial_x[offset + 2] = extrinsic_params.rotation()[2];
      initial_x[offset + 3] = t[0];
      initial_x[offset + 4] = t[1];
      initial_x[offset + 5] = t[2];
    }
    Index intrinsic_begin = vector_function.GetPointParamsSize() +
                            vector_function.GetExtrinsicParamsSize();
    Index intrinsic_params_size_per_camera =
      vector_function.GetIntrinsicParamsSizePerCamera();
    for (Index camera_id = 0; camera_id < number_of_cameras; camera_id++)
    {
      Index offset = intrinsic_begin +
                     camera_id * intrinsic_params_size_per_camera;
      const IntrinsicParams& intrinsic_params = intrinsic_params_set[camera_id];
      initial_x[offset + 0] = intrinsic_params.k1();
      initial_x[offset + 1] = intrinsic_params.k2();
      initial_x[offset + 2] = intrinsic_params.k3();
      initial_x[offset + 3] = intrinsic_params.d1();
      initial_x[offset + 4] = intrinsic_params.d2();
      initial_x[offset + 5] = intrinsic_params.focal_length();
      initial_x[offset + 6] = intrinsic_params.skew();
      initial_x[offset + 7] = intrinsic_params.principal_point_x();
      initial_x[offset + 8] = intrinsic_params.principal_point_y();
      initial_x[offset + 9] = intrinsic_params.pixel_ratio();
    }

    YVector near_y(y_size);
    Index feature_id = 0;
    for (size_t track_bundle_id = 0; track_bundle_id < number_of_tracks_bundle;
         track_bundle_id++)
    {
      size_t number_of_views_bundle = tracks_bundle[track_bundle_id].size();
      for (size_t view_id = 0; view_id < number_of_views_bundle; view_id++)
      {
        size_t image_id = tracks_bundle[track_bundle_id][view_id].first;
        size_t key_id = tracks_bundle[track_bundle_id][view_id].second;

        near_y.segment(feature_id * BAVectorFunction::params_per_key_,
                       BAVectorFunction::params_per_key_) =
          image_keysets[image_id][key_id];

        feature_id++;
      }
    }
    for (size_t gcp_estimate_id = 0; gcp_estimate_id < number_of_gcps_estimate;
         gcp_estimate_id++)
    {
      size_t number_of_views_gcp = tracks_gcp_estimate[gcp_estimate_id].size();
      for (size_t view_id = 0; view_id < number_of_views_gcp; view_id++)
      {
        size_t image_id = tracks_gcp_estimate[gcp_estimate_id][view_id].first;
        size_t key_id = tracks_gcp_estimate[gcp_estimate_id][view_id].second;

        near_y.segment(feature_id * BAVectorFunction::params_per_key_,
                       BAVectorFunction::params_per_key_) =
          image_keysets_gcp[image_id][key_id];
        feature_id++;
      }
    }
    Index y_offset = vector_function.GetYKeysSize();
    for (size_t point_constraint_id = 0;
         point_constraint_id < vector_function.point_constraints().size();
         point_constraint_id++)
    {
      const hs::sfm::ba::PointConstraint& point_constraint =
        vector_function.point_constraints()[point_constraint_id];
      size_t gcp_estimate_id = point_constraint.point_id -
                              number_of_tracks_bundle;
      size_t gcp_measure_id = estimate_measure_map[gcp_estimate_id];
      if (point_constraint.mask[hs::sfm::ba::POINT_CONSTRAIN_X])
      {
        near_y[y_offset] = gcps_measure[gcp_measure_id][0];
        y_offset++;
      }
      if (point_constraint.mask[hs::sfm::ba::POINT_CONSTRAIN_Y])
      {
        near_y[y_offset] = gcps_measure[gcp_measure_id][1];
        y_offset++;
      }
      if (point_constraint.mask[hs::sfm::ba::POINT_CONSTRAIN_Z])
      {
        near_y[y_offset] = gcps_measure[gcp_measure_id][2];
        y_offset++;
      }
    }

    YCovarianceInverse y_covariance_inverse;
    y_covariance_inverse.SetKeysUniformStdDev(image_accuracy_, number_of_keys);
    for (size_t i = 0; i < number_of_feature_map_gcp; i++)
    {
      size_t gcp_key_id = number_of_keys - number_of_feature_map_gcp + i;
      KeyBlock& key_block = y_covariance_inverse.GetKeyBlock(gcp_key_id);
      key_block(0, 0) = Scalar(gcp_image_accuracy_);
      key_block(1, 1) = Scalar(gcp_image_accuracy_);
    }
    auto itr_point_constriant = vector_function.point_constraints().begin();
    auto itr_point_constriant_end = vector_function.point_constraints().end();
    for (; itr_point_constriant != itr_point_constriant_end;
         ++itr_point_constriant)
    {
      Scalar planar_constraint = Scalar(1) / (gcp_planar_accuracy_);
      Scalar height_constraint = Scalar(1) / (gcp_height_accuracy_);
      if (itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_X] &&
          itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Y] &&
          !itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Z])
      {
        y_covariance_inverse.AddConstraint(planar_constraint);
        y_covariance_inverse.AddConstraint(planar_constraint);
      }
      else if (itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_X] &&
               itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Y] &&
               itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Z])
      {
        y_covariance_inverse.AddConstraint(planar_constraint);
        y_covariance_inverse.AddConstraint(planar_constraint);
        y_covariance_inverse.AddConstraint(height_constraint);
      }
    }

    int max_num_iterations = 50;
    double function_tolerance = 1e-8;
    double parameter_tolerance = 1e-10;
    BAOptimizor ba_optimizor(initial_x, int(number_of_threads_),
                             max_num_iterations, function_tolerance,
                             parameter_tolerance);
    XVector optimized_x;
    if (ba_optimizor(vector_function, near_y, y_covariance_inverse,
                     optimized_x) != 0)
    {
      return -1;
    }

    for (size_t track_bundle_id = 0; track_bundle_id < number_of_tracks_bundle;
         track_bundle_id++)
    {
      size_t point_id = track_point_map_bundle[track_bundle_id];
      points[point_id] = vector_function.GetPoint(Index(track_bundle_id),
                                                  optimized_x);
    }
    for (size_t gcp_id = 0; gcp_id < number_of_gcps; gcp_id++)
    {
      Index x_id = Index(number_of_tracks_bundle + gcp_id);
      gcps_estimate[gcp_id] = vector_function.GetPoint(x_id, optimized_x);
    }
    for (Index image_id = 0; image_id < number_of_images; image_id++)
    {
      ExtrinsicParams& extrinsic_params = extrinsic_params_set[image_id];
      BAImage ba_image = vector_function.GetImage(image_id, optimized_x);
      extrinsic_params.rotation()[0] = ba_image[0];
      extrinsic_params.rotation()[1] = ba_image[1];
      extrinsic_params.rotation()[2] = ba_image[2];
      Position position = ba_image.template segment<3>(3);
      extrinsic_params.position() =
        -(extrinsic_params.rotation().Inverse() * position);
    }

    for (Index camera_id = 0; camera_id < number_of_cameras; camera_id++)
    {
      IntrinsicParams& intrinsic_params = intrinsic_params_set[camera_id];
      BACamera ba_camera = vector_function.GetCamera(camera_id, optimized_x);
      intrinsic_params.set_k1(ba_camera[0]);
      intrinsic_params.set_k2(ba_camera[1]);
      intrinsic_params.set_k3(ba_camera[2]);
      intrinsic_params.set_d1(ba_camera[3]);
      intrinsic_params.set_d2(ba_camera[4]);
      intrinsic_params.set_focal_length(ba_camera[5]);
      intrinsic_params.set_skew(ba_camera[6]);
      intrinsic_params.set_principal_point_x(ba_camera[7]);
      intrinsic_params.set_principal_point_y(ba_camera[8]);
      intrinsic_params.set_pixel_ratio(ba_camera[9]);
    }

    return 0;
  }

private:
  size_t number_of_threads_;
  Scalar gcp_planar_accuracy_;
  Scalar gcp_height_accuracy_;
  Scalar image_accuracy_;
  Scalar gcp_image_accuracy_;
};

}
}
}

#endif
