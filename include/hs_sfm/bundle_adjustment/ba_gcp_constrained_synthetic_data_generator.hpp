#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_SYNTHETIC_DATA_GENERATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_SYNTHETIC_DATA_GENERATOR_HPP_

#include "hs_sfm/utility/synthetic_scene_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar, typename _ImageDimension>
class BAGCPConstrainedSyntheticDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;

  typedef int Err;

  typedef hs::sfm::SceneGenerator<Scalar, ImageDimension> SceneGenerator;
  typedef typename SceneGenerator::IntrinsicParams IntrinsicParams;
  typedef typename SceneGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename SceneGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename SceneGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename SceneGenerator::Image Image;
  typedef typename SceneGenerator::ImageContainer ImageContainer;
  typedef typename SceneGenerator::Point3D Point3D;
  typedef typename SceneGenerator::Point3DContainer Point3DContainer;

  typedef hs::sfm::KeysGenerator<Scalar, ImageDimension> KeysGenerator;
  typedef typename KeysGenerator::Keys Keys;
  typedef typename KeysGenerator::KeysContainer KeysContainer;
  typedef typename KeysGenerator::Track Track;
  typedef typename KeysGenerator::TrackContainer TrackContainer;
  typedef typename KeysGenerator::CameraViewContainer CameraViewContainer;

  typedef BAGCPConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::FeatureMap FeatureMap;
  typedef typename VectorFunction::FeatureMapContainer FeatureMapContainer;

private:
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;

public:
  BAGCPConstrainedSyntheticDataGenerator(
    Scalar focal_length_in_metre,
    size_t number_of_strps,
    size_t number_of_cameras_in_strip,
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    size_t number_of_points,
    Scalar lateral_overlap_ratio,
    Scalar longitudinal_overlap_ratio,
    Scalar scene_max_height,
    Scalar camera_height_stddev,
    Scalar camera_plannar_stddev,
    Scalar camera_rot_stddev,
    Scalar north_west_angle,
    size_t number_of_gcps)
    : scene_generator_(focal_length_in_metre,
                       number_of_strps,
                       number_of_cameras_in_strip,
                       ground_resolution,
                       image_width,
                       image_height,
                       pixel_size,
                       number_of_points + number_of_gcps,
                       lateral_overlap_ratio,
                       longitudinal_overlap_ratio,
                       scene_max_height,
                       camera_height_stddev,
                       camera_plannar_stddev,
                       camera_rot_stddev,
                       north_west_angle),
      keys_generator_(image_width, image_height),
      number_of_gcps_(number_of_gcps) {}

  Err operator() (VectorFunction& vector_function,
                  XVector& x, YVector& y) const
  {
    //Copy 从ba_naive_synthetic_data_generator.hpp拷贝了大量代码过来
    //TODO:这里应该要有更合适的设计才对！
    IntrinsicParamsContainer intrinsic_params_set;
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    Point3DContainer points;
    if (scene_generator_(intrinsic_params_set,
                         extrinsic_params_set,
                         images,
                         points) != 0) 
      return -1;
    Scalar f = GetFocalLengthInPixel();

    KeysContainer keys_set;
    TrackContainer tracks;
    CameraViewContainer camera_views;
    if (keys_generator_(intrinsic_params_set,
                        extrinsic_params_set,
                        points, keys_set,
                        tracks,
                        camera_views) != 0) 
      return -1;

    size_t number_of_keys = 0;
    auto itr_keys = keys_set.begin();
    auto itr_keys_end = keys_set.end();
    std::vector<Index> keys_id_offsets(extrinsic_params_set.size());
    Index i = 0;
    for (; itr_keys != itr_keys_end; ++itr_keys, ++i)
    {
      keys_id_offsets[i] = number_of_keys;
      number_of_keys += itr_keys->size();
    }

    std::vector<size_t> camera_map(extrinsic_params_set.size(), 0);
    auto itr_camera_view = camera_views.begin();
    auto itr_camera_view_end = camera_views.end();
    i = 0;
    size_t number_of_view_key_cameras = 0;
    for (; itr_camera_view != itr_camera_view_end; ++itr_camera_view, ++i)
    {
      if (!itr_camera_view->empty())
      {
        camera_map[i] = number_of_view_key_cameras + 1;
        number_of_view_key_cameras++;
      }
    }

    std::vector<size_t> point_map(points.size(), 0);
    auto itr_track = tracks.begin();
    auto itr_track_end = tracks.end();
    i = 0;
    size_t number_of_view_key_point = 0;
    for (; itr_track != itr_track_end; ++itr_track, ++i)
    {
      if (!itr_track->empty())
      {
        point_map[i] = number_of_view_key_point + 1;
        number_of_view_key_point++;
      }
    }

    itr_track = tracks.begin();
    itr_track_end = tracks.end();
    FeatureMapContainer feature_maps(number_of_keys);
    i = 0;
    for (; itr_track != itr_track_end; ++itr_track, ++i)
    {
      auto itr_view = itr_track->begin();
      auto itr_view_end = itr_track->end();
      for (; itr_view != itr_view_end; ++itr_view)
      {
        Index feature_id = keys_id_offsets[itr_view->first] + 
                           itr_view->second;
        feature_maps[feature_id].first = Index(camera_map[itr_view->first] - 1);
        feature_maps[feature_id].second = Index(point_map[i] - 1);
      }
    }

    Index number_of_cameras = Index(number_of_view_key_cameras);
    Index number_of_points = Index(number_of_view_key_point);
    Index number_of_features = Index(number_of_keys);
    vector_function.set_number_of_cameras(number_of_cameras);
    vector_function.set_number_of_points(number_of_points);
    vector_function.set_number_of_features(number_of_features);
    vector_function.set_feature_maps(feature_maps);
    vector_function.set_number_of_gcps(Index(number_of_gcps_));

    //相机外参数
    Index x_size = vector_function.GetXSize();
    x.resize(x_size);
    Index y_size = vector_function.GetYSize();
    y.resize(y_size);

    auto ext_itr = extrinsic_params_set.begin();
    auto ext_itr_end = extrinsic_params_set.end();
    itr_camera_view = camera_views.begin();
    itr_camera_view_end = camera_views.end();
    i = 0;
    for (; ext_itr != ext_itr_end; ++ext_itr, ++itr_camera_view)
    {
      if (!itr_camera_view->empty())
      {
        Vector3 t = -(ext_itr->rotation() * ext_itr->position());
        x[i * VectorFunction::params_per_camera_ + 0] =
          ext_itr->rotation()[0];
        x[i * VectorFunction::params_per_camera_ + 1] =
          ext_itr->rotation()[1];
        x[i * VectorFunction::params_per_camera_ + 2] =
          ext_itr->rotation()[2];
        x[i * VectorFunction::params_per_camera_ + 3] = t[0];
        x[i * VectorFunction::params_per_camera_ + 4] = t[1];
        x[i * VectorFunction::params_per_camera_ + 5] = t[2];
        i++;
      }
    }

    //点云
    Index camera_params_size = vector_function.GetCameraParamsSize();
    auto point_itr = points.begin();
    auto point_itr_end = points.end();
    itr_track = tracks.begin();
    itr_track_end = tracks.end();
    i = 0;
    for (; point_itr != point_itr_end; ++point_itr, itr_track++)
    {
      if (!itr_track->empty())
      {
        x.segment(camera_params_size + i * VectorFunction::params_per_point_,
                  VectorFunction::params_per_point_) = *point_itr;
        i++;
      }
    }

    //特征点
    itr_keys = keys_set.begin();
    i = 0;
    for (; itr_keys != itr_keys_end; ++itr_keys)
    {
      for (size_t keyId = 0; keyId < itr_keys->size(); keyId++)
      {
        y.segment(i * VectorFunction::params_per_feature_,
                  VectorFunction::params_per_feature_) =
          (*itr_keys)[keyId] / f;
        i++;
      }
    }

    //像控点
    Index y_gcp_start = y_size - Index(number_of_gcps_) *
                                 VectorFunction::params_per_point_;
    Index x_gcp_start = x_size - Index(number_of_gcps_) *
                                 VectorFunction::params_per_point_;

    for (Index i = 0; i < Index(number_of_gcps_); i++)
    {
      y.segment(y_gcp_start + i * VectorFunction::params_per_point_,
        VectorFunction::params_per_point_) =
        x.segment(x_gcp_start + i * VectorFunction::params_per_point_,
        VectorFunction::params_per_point_);
    }

    return 0;
  }

  inline Scalar GetFocalLengthInPixel() const
  {
    return scene_generator_.GetFocalLengthInPixel();
  }

private:
  SceneGenerator scene_generator_;
  KeysGenerator keys_generator_;
  size_t number_of_gcps_;
};

}
}
}

#endif