#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_SYNTHETIC_DATA_GENERATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_SYNTHETIC_DATA_GENERATOR_HPP_

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/sfm_utility/synthetic_scene_generator.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar, typename _ImageDimension>
class BANaiveSyntheticDataGenerator
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

  typedef BANaiveVectorFunction<Scalar> BAVectorFunction;
  typedef typename BAVectorFunction::Index Index;
  typedef typename BAVectorFunction::XVector XVector;
  typedef typename BAVectorFunction::YVector YVector;
  typedef typename BAVectorFunction::FeatureMap FeatureMap;
  typedef typename BAVectorFunction::FeatureMapContainer FeatureMapContainer;
  typedef typename BAVectorFunction::Vector3 Vector3;

  BANaiveSyntheticDataGenerator(
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
    Scalar north_west_angle)
    : scene_generator_(focal_length_in_metre,
                       number_of_strps,
                       number_of_cameras_in_strip,
                       ground_resolution,
                       image_width,
                       image_height,
                       pixel_size,
                       number_of_points,
                       lateral_overlap_ratio,
                       longitudinal_overlap_ratio,
                       scene_max_height,
                       camera_height_stddev,
                       camera_plannar_stddev,
                       camera_rot_stddev,
                       north_west_angle),
      keys_generator_(image_width, image_height) {}

  Err operator () (BAVectorFunction& vector_function,
                   XVector& x, YVector& y) const
  {
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
        x[i * BAVectorFunction::params_per_camera_ + 0] =
          ext_itr->rotation()[0];
        x[i * BAVectorFunction::params_per_camera_ + 1] =
          ext_itr->rotation()[1];
        x[i * BAVectorFunction::params_per_camera_ + 2] =
          ext_itr->rotation()[2];
        x[i * BAVectorFunction::params_per_camera_ + 3] = t[0];
        x[i * BAVectorFunction::params_per_camera_ + 4] = t[1];
        x[i * BAVectorFunction::params_per_camera_ + 5] = t[2];
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
        x.segment(camera_params_size + i * BAVectorFunction::params_per_point_,
                  BAVectorFunction::params_per_point_) = *point_itr;
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
        y.segment(i * BAVectorFunction::params_per_feature_,
                  BAVectorFunction::params_per_feature_) =
          (*itr_keys)[keyId] / f;
        i++;
      }
    }

    return 0;
  }

  inline Scalar GetFocalLengthInPixel() const
  {
    return scene_generator_.GetFocalLengthInPixel();
  }

  Err GenerateGCPs(size_t number_of_gcps,
                   const BAVectorFunction& vector_function,
                   const XVector& x,
                   Point3DContainer& gcps,
                   KeysContainer& image_keys_set,
                   TrackContainer& tracks) const
  {
    scene_generator_.GenerateScenePoints(number_of_gcps, gcps);
    Index number_of_cameras = vector_function.number_of_cameras();

    IntrinsicParamsContainer intrinsic_params_set(
      number_of_cameras, IntrinsicParams(GetFocalLengthInPixel()));
    ExtrinsicParamsContainer extrinsic_params_set(number_of_cameras);

    for (Index i = 0; i < number_of_cameras; i++)
    {
      extrinsic_params_set[i].rotation()[0] =
        x[i * BAVectorFunction::params_per_camera_ + 0];
      extrinsic_params_set[i].rotation()[1] =
        x[i * BAVectorFunction::params_per_camera_ + 1];
      extrinsic_params_set[i].rotation()[2] =
        x[i * BAVectorFunction::params_per_camera_ + 2];
      Vector3 t = x.segment(i * BAVectorFunction::params_per_camera_ + 3, 3);
      extrinsic_params_set[i].position() =
        -(extrinsic_params_set[i].rotation().Inverse() * t);
    }

    CameraViewContainer camera_views;
    return keys_generator_(intrinsic_params_set,
                           extrinsic_params_set,
                           gcps,
                           image_keys_set,
                           tracks,
                           camera_views);
  }

private:
  SceneGenerator scene_generator_;
  KeysGenerator keys_generator_;
};

}
}
}

#endif
