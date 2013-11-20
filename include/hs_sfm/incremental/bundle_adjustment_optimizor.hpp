#ifndef _HS_SFM_INCREMENTAL_BUNDLE_ADJUSTMENT_OPTIMIZOR_HPP_
#define _HS_SFM_INCREMENTAL_BUNDLE_ADJUSTMENT_OPTIMIZOR_HPP_

#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_levenberg_marquardt_optimizor.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar>
class BundleAdjustmentOptimizor
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
  typedef ObjectIndexMap ImageExtrinsicMap;

private:
  typedef hs::sfm::ba::BANaiveVectorFunction<Scalar> BAVectorFunction;
  typedef typename BAVectorFunction::Index Index;
  typedef typename BAVectorFunction::XVector XVector;
  typedef typename BAVectorFunction::YVector YVector;
  typedef typename BAVectorFunction::FeatureMap FeatureMap;
  typedef typename BAVectorFunction::FeatureMapContainer FeatureMapContainer;

  typedef hs::sfm::ba::BANaiveLevenbergMarquardtOptimizor<BAVectorFunction>
          BAOptimizor;
  typedef typename BAOptimizor::YCovarianceInverse YCovarianceInverse;

  typedef typename IntrinsicParams::KMatrix KMatrix;

public:
  Err operator() (const ImageKeysetContainer& image_keysets,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const TrackContainer& tracks,
                  const ImageExtrinsicMap& image_extrinsic_map,
                  const TrackPointMap& track_point_map,
                  const ViewInfoIndexer& view_info_indexer,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  PointContainer& points) const
  {
    size_t number_of_tracks = tracks.size();
    TrackContainer tracks_bundle;
    ObjectIndexMap bundle_point_map;
    Index number_of_features = 0;
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
        }
        if (!track_bundle.empty())
        {
          tracks_bundle.push_back(track_bundle);
          size_t point_id = track_point_map[track_id];
          bundle_point_map.AddObject(point_id);
          number_of_features += Index(track_bundle.size());
        }
      }
    }

    size_t number_of_tracks_bundle = tracks_bundle.size();
    FeatureMapContainer feature_maps;
    for (size_t i = 0; i < number_of_tracks_bundle; i++)
    {
      size_t number_of_views_bundle = tracks_bundle[i].size();
      for (size_t j = 0; j < number_of_views_bundle; j++)
      {
        size_t image_id = tracks_bundle[i][j].first;
        size_t key_id = tracks_bundle[i][j].second;
        size_t extrinsic_id = image_extrinsic_map[image_id];

        feature_maps.push_back(std::make_pair(extrinsic_id, i));
      }
    }

    Index number_of_points = Index(number_of_tracks_bundle);
    Index number_of_cameras = Index(extrinsic_params_set.size());
    BAVectorFunction vector_function;
    vector_function.set_number_of_cameras(number_of_cameras);
    vector_function.set_number_of_points(number_of_points);
    vector_function.set_number_of_features(number_of_features);
    vector_function.set_feature_maps(feature_maps);

    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    XVector initial_x(x_size);
    YVector near_y(y_size);
    for (Index i = 0; i < number_of_cameras; i++)
    {
      const ExtrinsicParams& extrinsic_params = extrinsic_params_set[i];
      Point t = -(extrinsic_params.rotation() * extrinsic_params.position());
      initial_x[i * BAVectorFunction::params_per_camera_ + 0] =
        extrinsic_params.rotation()[0];
      initial_x[i * BAVectorFunction::params_per_camera_ + 1] =
        extrinsic_params.rotation()[1];
      initial_x[i * BAVectorFunction::params_per_camera_ + 2] =
        extrinsic_params.rotation()[2];
      initial_x[i * BAVectorFunction::params_per_camera_ + 3] = t[0];
      initial_x[i * BAVectorFunction::params_per_camera_ + 4] = t[1];
      initial_x[i * BAVectorFunction::params_per_camera_ + 5] = t[2];
    }
    Index camera_params_size = vector_function.GetCameraParamsSize();
    for (Index i = 0; i < number_of_points; i++)
    {
      size_t point_id = bundle_point_map[size_t(i)];
      initial_x.segment(
        camera_params_size + i * BAVectorFunction::params_per_point_,
        BAVectorFunction::params_per_point_) = points[point_id];
    }
    Index feature_id = 0;
    for (size_t i = 0; i < number_of_tracks_bundle; i++)
    {
      size_t number_of_views_bundle = tracks_bundle[i].size();
      for (size_t j = 0; j < number_of_views_bundle; j++)
      {
        size_t image_id = tracks_bundle[i][j].first;
        size_t key_id = tracks_bundle[i][j].second;
        size_t extrinsic_id = image_extrinsic_map[image_id];

        KMatrix K = intrinsic_params_set[image_id].GetKMatrix();
        EIGEN_VECTOR(Scalar, 3) hkey;
        hkey.template segment<2>(0) = image_keysets[image_id][key_id];
        hkey(2) = Scalar(1);
        hkey = K.inverse() * hkey;
        hkey /= hkey(2);

        near_y.segment(feature_id * BAVectorFunction::params_per_feature_,
                       BAVectorFunction::params_per_feature_) =
          hkey.segment(0, 2);
        feature_id++;
      }
    }
    YCovarianceInverse y_covariance_inverse;
    for (Index i = 0; i < number_of_features; i++)
    {
      y_covariance_inverse.blocks.push_back(
        EIGEN_MATRIX(Scalar, 2, 2)::Identity());
    }

    BAOptimizor ba_optimizor(initial_x);
    XVector optimized_x;
    if (ba_optimizor(vector_function, near_y, y_covariance_inverse,
                     optimized_x) != 0)
    {
      return -1;
    }

    for (Index i = 0; i < number_of_cameras; i++)
    {
      ExtrinsicParams& extrinsic_params = extrinsic_params_set[i];
      extrinsic_params.rotation()[0] =
        optimized_x[i * BAVectorFunction::params_per_camera_ + 0];
      extrinsic_params.rotation()[1] =
        optimized_x[i * BAVectorFunction::params_per_camera_ + 1];
      extrinsic_params.rotation()[2] =
        optimized_x[i * BAVectorFunction::params_per_camera_ + 2];
      Point t =
        optimized_x.segment(i * BAVectorFunction::params_per_camera_ + 3, 3);
      Point c = -((extrinsic_params.rotation().Inverse()) * t);
      extrinsic_params.position() = c;
    }

    for (Index i = 0; i < number_of_points; i++)
    {
      size_t point_id = bundle_point_map[size_t(i)];
      points[point_id] =
        optimized_x.segment(
        camera_params_size + i * BAVectorFunction::params_per_point_,
        BAVectorFunction::params_per_point_);
    }

    return 0;
  }
};

}
}
}

#endif
