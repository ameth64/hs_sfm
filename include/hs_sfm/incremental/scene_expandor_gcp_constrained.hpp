#ifndef _HS_SFM_INCREMENTAL_SCENE_EXPANDOR_GCP_CONSTRAINED_HPP_
#define _HS_SFM_INCREMENTAL_SCENE_EXPANDOR_GCP_CONSTRAINED_HPP_

#include "hs_sfm/sfm_utility/similar_transform_estimator.hpp"
#include "hs_sfm/incremental/scene_expandor.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_levenberg_marquardt_optimizor.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar>
class SceneExpandorGCPConstrained : public SceneExpandor<_Scalar>
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

private:
  typedef SceneExpandor<Scalar> Base;
  using Base::Initialize;
  using Base::AddNewImage;
  using Base::TriangulateNewPoints;
  using Base::BundleAdjustment;
  using Base::GetReprojectiveError;
  using Base::min_triangulate_views_;
  using Base::triangulate_error_threshold_;

  typedef hs::sfm::SimilarTransformEstimator<Scalar>
          TransformEstimator;
  typedef typename TransformEstimator::Rotation Rotation;
  typedef typename TransformEstimator::Translate Translate;

  typedef hs::sfm::ba::BAGCPConstrainedVectorFunction<Scalar> GCPVectorFunction;
  typedef typename GCPVectorFunction::Index Index;
  typedef typename GCPVectorFunction::XVector GCPXVector;
  typedef typename GCPVectorFunction::YVector GCPYVector;
  typedef typename GCPVectorFunction::FeatureMap FeatureMap;
  typedef typename GCPVectorFunction::FeatureMapContainer FeatureMapContainer;

  typedef hs::sfm::ba::BAGCPConstrainedLevenbergMarquardtOptimizor<
            GCPVectorFunction> GCPOptimizor;
  typedef typename GCPOptimizor::YCovarianceInverse GCPYCovarianceInverse;

  typedef typename Base::ViewInfoIndexer ViewInfoIndexer;
  typedef typename Base::ImageViewTracksContainer
          ImageViewTracksContainer;
  typedef typename Base::ExtrinsicParams ExtrinsicParams;
  typedef typename Base::Point Point;
  typedef typename Base::ViewInfo ViewInfo;
  typedef typename Base::KMatrix KMatrix;

public:
  typedef typename Base::ImageKeysetContainer ImageKeysetContainer;
  typedef typename Base::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef typename Base::PointContainer PointContainer;
  typedef typename Base::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef typename Base::ImageExtrinsicMap ImageExtrinsicMap;
  typedef typename Base::TrackPointMap TrackPointMap;

public:
  SceneExpandorGCPConstrained(
    size_t add_new_image_matches_threshold = 8,
    Scalar pmatrix_ransac_threshold = 4.0,
    size_t min_triangulate_views = 2,
    Scalar triangulate_error_threshold = 4.0)
    : Base(add_new_image_matches_threshold,
           pmatrix_ransac_threshold,
           min_triangulate_views,
           triangulate_error_threshold) {}

public:
  Err operator() (const ImageKeysetContainer& image_keysets,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const hs::sfm::TrackContainer& tracks,
                  const PointContainer& gcps,
                  const ImageKeysetContainer& image_key_sets_gcp,
                  const hs::sfm::TrackContainer& tracks_gcp,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  ImageExtrinsicMap& image_extrinsic_map,
                  PointContainer& points,
                  TrackPointMap& track_point_map,
                  Scalar& reprojection_error) const
  {
    return Run(image_keysets,
               intrinsic_params_set,
               tracks,
               gcps,
               image_key_sets_gcp,
               tracks_gcp,
               extrinsic_params_set,
               image_extrinsic_map,
               points,
               track_point_map,
               reprojection_error);
  }

private:
  Err Run(const ImageKeysetContainer& image_keysets,
          const IntrinsicParamsContainer& intrinsic_params_set,
          const hs::sfm::TrackContainer& tracks,
          const PointContainer& gcps,
          const ImageKeysetContainer& image_key_sets_gcp,
          const hs::sfm::TrackContainer& tracks_gcp,
          ExtrinsicParamsContainer& extrinsic_params_set,
          ImageExtrinsicMap& image_extrinsic_map,
          PointContainer& points,
          TrackPointMap& track_point_map,
          Scalar& reprojection_error) const
  {
    ImageViewTracksContainer image_view_tracks_set;
    ViewInfoIndexer view_info_indexer;
    if (Initialize(image_keysets,
      intrinsic_params_set,
      tracks,
      image_extrinsic_map,
      track_point_map,
      image_view_tracks_set,
      view_info_indexer) != 0)
    {
      return -1;
    }

    while (1)
    {
      Scalar reprojection_error_before =
      GetReprojectiveError(image_keysets,
                           intrinsic_params_set,
                           tracks,
                           image_extrinsic_map,
                           track_point_map,
                           view_info_indexer,
                           extrinsic_params_set,
                           points);
      reprojection_error = reprojection_error_before;
      if (BundleAdjustment(image_keysets,
                           intrinsic_params_set,
                           tracks,
                           image_extrinsic_map,
                           track_point_map,
                           view_info_indexer,
                           extrinsic_params_set,
                           points) != 0)
      {
        break;
      }
      reprojection_error =
        GetReprojectiveError(image_keysets,
                             intrinsic_params_set,
                             tracks,
                             image_extrinsic_map,
                             track_point_map,
                             view_info_indexer,
                             extrinsic_params_set,
                             points);

      ExtrinsicParams new_extrinsic_params;
      size_t new_image_id;
      if (AddNewImage(image_keysets,
                      intrinsic_params_set,
                      points,
                      track_point_map,
                      image_view_tracks_set,
                      image_extrinsic_map,
                      view_info_indexer,
                      new_extrinsic_params,
                      new_image_id) != 0)
      {
        break;
      }

      extrinsic_params_set.push_back(new_extrinsic_params);
      image_extrinsic_map[new_image_id] = extrinsic_params_set.size() - 1;

      if (TriangulateNewPoints(image_keysets,
                               intrinsic_params_set,
                               tracks,
                               extrinsic_params_set,
                               image_extrinsic_map,
                               min_triangulate_views_,
                               triangulate_error_threshold_,
                               points,
                               track_point_map,
                               view_info_indexer) != 0)
      {
        break;
      }
    }

    if (BundleAdjustmentWithGCPConstraits(image_keysets,
                                          intrinsic_params_set,
                                          tracks,
                                          image_extrinsic_map,
                                          track_point_map,
                                          view_info_indexer,
                                          gcps,
                                          image_key_sets_gcp,
                                          tracks_gcp,
                                          extrinsic_params_set,
                                          points) != 0)
    {
      return -1;
    }

    reprojection_error = GetReprojectiveError(image_keysets,
                                              intrinsic_params_set,
                                              tracks,
                                              image_extrinsic_map,
                                              track_point_map,
                                              view_info_indexer,
                                              extrinsic_params_set,
                                              points);

    return 0;
  }

private:
  Err BundleAdjustmentWithGCPConstraits(
    const ImageKeysetContainer& image_keysets,
    const IntrinsicParamsContainer& intrinsic_params_set,
    const hs::sfm::TrackContainer& tracks,
    const ImageExtrinsicMap& image_extrinsic_map,
    const TrackPointMap& track_point_map,
    const ViewInfoIndexer& view_info_indexer,
    const PointContainer& gcps,
    const ImageKeysetContainer& image_key_sets_gcp,
    const hs::sfm::TrackContainer& tracks_gcp,
    ExtrinsicParamsContainer& extrinsic_params_set,
    PointContainer& points) const
  {
    //计算像控点在相对坐标系下的坐标
    size_t number_of_gcps = gcps.size();
    TrackPointMap track_point_map_gcp(number_of_gcps);
    ViewInfoIndexer gcp_view_info_indexer;
    gcp_view_info_indexer.SetViewInfoByTracks(tracks_gcp,
                                              track_point_map_gcp,
                                              image_extrinsic_map);

    PointContainer gcps_rel;
    if (TriangulateNewPoints(image_key_sets_gcp,
                             intrinsic_params_set,
                             tracks_gcp,
                             extrinsic_params_set,
                             image_extrinsic_map,
                             min_triangulate_views_,
                             triangulate_error_threshold_ * Scalar(2),
                             gcps_rel,
                             track_point_map_gcp,
                             gcp_view_info_indexer) != 0)
    {
      return -1;
    }

    //计算相似变换，将相对坐标系转到绝对坐标系
    size_t number_of_available_gcps = gcps_rel.size();
    if (number_of_available_gcps < 4)
    {
      return -1;
    }
    PointContainer gcps_abs;
    for (size_t i = 0 ; i < number_of_gcps; i++)
    {
      if (track_point_map_gcp.IsValid(i))
      {
        gcps_abs.push_back(gcps[i]);
      }
    }

    TransformEstimator transform_estimator;
    Rotation rotation_similar;
    Translate translate_similar;
    Scalar scale_similar;
    if (transform_estimator(gcps_rel, gcps_abs,
                            rotation_similar,
                            translate_similar,
                            scale_similar) != 0)
    {
      return -1;
    }

    auto itr_extrinsic = extrinsic_params_set.begin();
    auto itr_extrinsic_end = extrinsic_params_set.end();
    for (; itr_extrinsic != itr_extrinsic_end; ++itr_extrinsic)
    {
      itr_extrinsic->rotation() =
        itr_extrinsic->rotation() * rotation_similar.Inverse();
      itr_extrinsic->position() =
        scale_similar * (rotation_similar * itr_extrinsic->position()) +
        translate_similar;
    }
    auto itr_point = points.begin();
    auto itr_point_end = points.end();
    for (; itr_point != itr_point_end; ++itr_point)
    {
      Point position_rel = *itr_point;
      (*itr_point) =
        scale_similar * (rotation_similar * position_rel) + translate_similar;
    }

    Scalar mean_gcp_error = Scalar(0);
    for (size_t i = 0; i < number_of_available_gcps; i++)
    {
      Point gcp_abs = gcps_rel[i];
      gcp_abs = scale_similar * (rotation_similar * gcp_abs) +
                translate_similar;
      Point diff = gcps_abs[i] - gcp_abs;
      Scalar error = diff.norm();
      mean_gcp_error += error;
    }
    mean_gcp_error /= Scalar(number_of_available_gcps);

    size_t number_of_tracks_gcp = tracks_gcp.size();
    TrackContainer tracks_bundle_gcp;
    ObjectIndexMap bundle_point_map_gcp;
    for (size_t track_id = 0; track_id < number_of_tracks_gcp; track_id++)
    {
      if (track_point_map_gcp.IsValid(track_id))
      {
        size_t number_of_views = tracks_gcp[track_id].size();
        hs::sfm::Track track_bundle_gcp;
        for (size_t i = 0; i < number_of_views; i++)
        {
          size_t image_id = tracks_gcp[track_id][i].first;
          const ViewInfo& view_info =
            gcp_view_info_indexer.GetViewInfoByTrackImage(track_id, image_id);
          if (image_extrinsic_map.IsValid(image_id) &&
              !view_info.is_blunder)
          {
            track_bundle_gcp.push_back(tracks_gcp[track_id][i]);
          }
        }//for (size_t i = 0; i < number_of_views; i++)

        if (!track_bundle_gcp.empty())
        {
          tracks_bundle_gcp.push_back(track_bundle_gcp);
          size_t point_id = track_point_map_gcp[track_id];
          bundle_point_map_gcp.AddObject(point_id);
        }//if (!track_bundle_gcp.empty())
      }//if (track_point_map_gcp.IsValid(track_id))
    }//for (size_t track_id = 0; track_id < number_of_tracks_gcp; track_id++)

    size_t number_of_tracks = tracks.size();
    TrackContainer tracks_bundle;
    ObjectIndexMap bundle_point_map;
    for (size_t track_id = 0; track_id < number_of_tracks; track_id++)
    {
      if (track_point_map.IsValid(track_id))
      {
        size_t number_of_views = tracks[track_id].size();
        hs::sfm::Track track_bundle;
        for (size_t i = 0; i < number_of_views; i++)
        {
          size_t image_id = tracks[track_id][i].first;
          const ViewInfo& view_info =
            view_info_indexer.GetViewInfoByTrackImage(track_id, image_id);
          if (image_extrinsic_map.IsValid(image_id) &&
              !view_info.is_blunder)
          {
            track_bundle.push_back(tracks[track_id][i]);
          }
        }//for (size_t i = 0; i < number_of_views; i++)

        if (!track_bundle.empty())
        {
          tracks_bundle.push_back(track_bundle);
          size_t point_id = track_point_map[track_id];
          bundle_point_map.AddObject(point_id);
        }
      }//if (track_point_map.IsValid(track_id))
    }//for (size_t track_id = 0; track_id < number_of_tracks; track_id++)

    FeatureMapContainer feature_maps;
    size_t number_of_tracks_bundle = tracks_bundle.size();
    for (size_t i = 0; i < number_of_tracks_bundle; i++)
    {
      size_t number_of_view_bundle = tracks_bundle[i].size();
      for (size_t j = 0; j < number_of_view_bundle; j++)
      {
        size_t image_id = tracks_bundle[i][j].first;
        size_t extrinsic_id = image_extrinsic_map[image_id];

        feature_maps.push_back(std::make_pair(extrinsic_id, i));
      }
    }

    size_t number_of_tracks_bundle_gcp = tracks_bundle_gcp.size();
    for (size_t i = 0; i < number_of_tracks_bundle_gcp; i++)
    {
      size_t number_of_views_bundle_gcp = tracks_bundle_gcp[i].size();
      for (size_t j = 0; j < number_of_views_bundle_gcp; j++)
      {
        size_t image_id = tracks_bundle_gcp[i][j].first;
        size_t extrinsic_id = image_extrinsic_map[image_id];

        feature_maps.push_back(std::make_pair(extrinsic_id,
                                              number_of_tracks_bundle + i));
      }
    }

    Index number_of_points = Index(number_of_tracks_bundle +
                                   number_of_tracks_bundle_gcp);
    Index number_of_cameras = Index(extrinsic_params_set.size());
    Index number_of_features = Index(feature_maps.size());
    GCPVectorFunction vector_function;
    vector_function.set_number_of_cameras(number_of_cameras);
    vector_function.set_number_of_points(number_of_points);
    vector_function.set_number_of_features(number_of_features);
    vector_function.set_number_of_gcps(Index(number_of_available_gcps));
    vector_function.set_feature_maps(feature_maps);

    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    GCPXVector initial_x(x_size);
    GCPYVector near_y(y_size);
    for (Index i = 0; i < number_of_cameras; i++)
    {
      const ExtrinsicParams& extrinsic_params = extrinsic_params_set[i];
      Point t = -(extrinsic_params.rotation() * extrinsic_params.position());
      initial_x[i * GCPVectorFunction::params_per_camera_ + 0] =
        extrinsic_params.rotation()[0];
      initial_x[i * GCPVectorFunction::params_per_camera_ + 1] =
        extrinsic_params.rotation()[1];
      initial_x[i * GCPVectorFunction::params_per_camera_ + 2] =
        extrinsic_params.rotation()[2];
      initial_x[i * GCPVectorFunction::params_per_camera_ + 3] = t[0];
      initial_x[i * GCPVectorFunction::params_per_camera_ + 4] = t[1];
      initial_x[i * GCPVectorFunction::params_per_camera_ + 5] = t[2];
    }
    Index camera_params_size = vector_function.GetCameraParamsSize();
    for (Index i = 0; i < Index(number_of_tracks_bundle); i++)
    {
      size_t point_id = bundle_point_map[size_t(i)];
      initial_x.segment(
        camera_params_size + i * GCPVectorFunction::params_per_point_,
        GCPVectorFunction::params_per_point_) = points[point_id];
    }
    for (Index i = 0; i < Index(number_of_tracks_bundle_gcp); i++)
    {
      Point gcp = gcps_rel[i];
      gcp = scale_similar * (rotation_similar * gcp) + translate_similar;
      size_t gcp_id = bundle_point_map_gcp[size_t(i)];
      initial_x.segment(
        camera_params_size +
        number_of_tracks_bundle * GCPVectorFunction::params_per_point_ +
        i * GCPVectorFunction::params_per_point_,
        GCPVectorFunction::params_per_point_) = gcp;
    }

    Index feature_id = 0;
    for (size_t i = 0; i < number_of_tracks_bundle; i++)
    {
      size_t number_of_views_bundle = tracks_bundle[i].size();
      for (size_t j = 0; j < number_of_views_bundle; j++)
      {
        size_t image_id = tracks_bundle[i][j].first;
        size_t key_id = tracks_bundle[i][j].second;

        KMatrix K = intrinsic_params_set[image_id].GetKMatrix();
        EIGEN_VECTOR(Scalar, 3) hkey;
        hkey.template segment<2>(0) = image_keysets[image_id][key_id];
        hkey(2) = Scalar(1);
        hkey = K.inverse() * hkey;
        hkey /= hkey(2);

        near_y.segment(feature_id * GCPVectorFunction::params_per_feature_,
                       GCPVectorFunction::params_per_feature_) =
          hkey.segment(0, 2);
        feature_id++;
      }
    }
    for (size_t i = 0; i < number_of_tracks_bundle_gcp; i++)
    {
      size_t number_of_views_bundle = tracks_bundle_gcp[i].size();
      for (size_t j = 0; j < number_of_views_bundle; j++)
      {
        size_t image_id = tracks_bundle_gcp[i][j].first;
        size_t key_id = tracks_bundle_gcp[i][j].second;

        KMatrix K = intrinsic_params_set[image_id].GetKMatrix();
        EIGEN_VECTOR(Scalar, 3) hkey;
        hkey.template segment<2>(0) = image_key_sets_gcp[image_id][key_id];
        hkey(2) = Scalar(1);
        hkey = K.inverse() * hkey;
        hkey /= hkey(2);

        near_y.segment(feature_id * GCPVectorFunction::params_per_feature_,
                       GCPVectorFunction::params_per_feature_) =
          hkey.segment(0, 2);
        feature_id++;
      }
    }
    Index feature_params_size =
      number_of_features * GCPVectorFunction::params_per_feature_;
    for (size_t i = 0; i < number_of_available_gcps; i++)
    {
      near_y.segment(feature_params_size +
                     i * GCPVectorFunction::params_per_point_,
                     GCPVectorFunction::params_per_point_) = gcps_abs[i];
    }
    
    GCPYCovarianceInverse y_covariance_inverse;
    for (Index i = 0; i < number_of_features; i++)
    {
      y_covariance_inverse.naive_y_covariance_inverse.blocks.push_back(
        EIGEN_MATRIX(Scalar, 2, 2)::Identity());
    }
    EIGEN_MATRIX(Scalar, 3, 3) gcp_covariance;
    gcp_covariance.setIdentity();
    for (Index i = 0; i < Index(number_of_available_gcps); i++)
    {
      y_covariance_inverse.gcp_blocks.push_back(gcp_covariance);
    }

    GCPOptimizor gcp_optimizor(initial_x);
    GCPXVector optimized_x;
    if (gcp_optimizor(vector_function, near_y, y_covariance_inverse,
                      optimized_x) != 0)
    {
      return -1;
    }

    for (Index i = 0; i < number_of_cameras; i++)
    {
      ExtrinsicParams& extrinsic_params = extrinsic_params_set[i];
      extrinsic_params.rotation()[0] =
        optimized_x[i * GCPVectorFunction::params_per_camera_ + 0];
      extrinsic_params.rotation()[1] =
        optimized_x[i * GCPVectorFunction::params_per_camera_ + 1];
      extrinsic_params.rotation()[2] =
        optimized_x[i * GCPVectorFunction::params_per_camera_ + 2];
      Point t =
        optimized_x.segment(i * GCPVectorFunction::params_per_camera_ + 3, 3);
      Point c = -((extrinsic_params.rotation().Inverse()) * t);
      extrinsic_params.position() = c;
    }

    for (Index i = 0; i < Index(number_of_tracks_bundle); i++)
    {
      size_t point_id = bundle_point_map[size_t(i)];
      points[point_id] = optimized_x.segment(
        camera_params_size + i * GCPVectorFunction::params_per_point_,
        GCPVectorFunction::params_per_point_);
    }

    return 0;
  }
};

}
}
}

#endif
