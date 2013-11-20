#ifndef _HS_SFM_UNIT_TEST_INCREMENTAL_SYNTHETIC_DATA_GENERATOR_HPP_
#define _HS_SFM_UNIT_TEST_INCREMENTAL_SYNTHETIC_DATA_GENERATOR_HPP_

#include <string>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"
#include "hs_math/geometry/euler_angles.hpp"

#include "hs_sfm/synthetic/keyset_generator.hpp"
#include "hs_sfm/synthetic/relative_generator.hpp"
#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar, typename _ImageDimension>
class SyntheticDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
  typedef hs::sfm::synthetic::SceneGenerator<Scalar, ImageDimension>
          SceneGenerator;

  typedef hs::sfm::synthetic::KeysetGenerator<Scalar, ImageDimension>
          KeysetGenerator;

  typedef hs::sfm::synthetic::RelativeGenerator<Scalar> RelativeGenerator;

public:
  typedef typename SceneGenerator::IntrinsicParams IntrinsicParams;
  typedef typename IntrinsicParams::KMatrix KMatrix;
  typedef typename SceneGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename SceneGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename SceneGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename SceneGenerator::Image Image;
  typedef typename SceneGenerator::ImageContainer ImageContainer;
  typedef typename SceneGenerator::Point3D Point3D;
  typedef typename SceneGenerator::Point3DContainer Point3DContainer;

  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

  typedef typename RelativeGenerator::RMatrix RMatrix;
  typedef typename RelativeGenerator::Translate Translate;

public:
  SyntheticDataGenerator(
    Scalar focal_length_in_metre,
    size_t number_of_strips,
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
    Scalar camera_rotation_stddev,
    Scalar north_west_angle,
    Scalar outlier_ratio,
    Scalar key_stddev)
  : scene_generator_(focal_length_in_metre,
                     number_of_strips,
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
                     camera_rotation_stddev,
                     north_west_angle),
    keyset_generator_(image_width, image_height),
    outlier_ratio_(outlier_ratio),
    key_stddev_(key_stddev) {}

public:
  Err GenerateAbsoluteScene(
    IntrinsicParamsContainer& intrinsic_params_set,
    ExtrinsicParamsContainer& extrinsic_params_set_absolute,
    ImageContainer& images,
    Point3DContainer& points_absolute,
    KeysetContainer& keysets,
    hs::sfm::TrackContainer& tracks,
    hs::sfm::CameraViewContainer& camera_views) const
  {
    if (scene_generator_(intrinsic_params_set,
                         extrinsic_params_set_absolute,
                         images,
                         points_absolute) != 0)
    {
      return -1;
    }

    if (keyset_generator_(intrinsic_params_set,
                          extrinsic_params_set_absolute,
                          points_absolute,
                          keysets,
                          tracks,
                          camera_views) != 0)
    {
      return -1;
    }

    return 0;
  }

  Err GenerateGCPData(
    const IntrinsicParamsContainer& intrinsic_params_set,
    const ExtrinsicParamsContainer& extrinsic_params_set_absolute,
    size_t number_of_gcps,
    KeysetContainer& keysets_gcp,
    Point3DContainer& gcps,
    hs::sfm::TrackContainer& tracks_gcp) const
  {
    scene_generator_.GenerateScenePoints(number_of_gcps, gcps);

    CameraViewContainer camera_view_gcp;
    if (keyset_generator_(intrinsic_params_set,
                          extrinsic_params_set_absolute,
                          gcps,
                          keysets_gcp,
                          tracks_gcp,
                          camera_view_gcp) != 0)
    {
      return -1;
    }

    return 0;
  }

  Err GenerateRelativeScene(
    const hs::sfm::TrackContainer& tracks,
    const hs::sfm::CameraViewContainer& camera_views,
    const ExtrinsicParamsContainer& extrinsic_params_set_absolute,
    const Point3DContainer& points_absolute,
    RMatrix& rotation_similar,
    Translate& translate_similar,
    Scalar& scale_similar,
    size_t& camera_id_identity,
    size_t& camera_id_relative,
    ExtrinsicParamsContainer& extrinsic_params_set_relative,
    Point3DContainer& points_relative) const
  {
    extrinsic_params_set_relative = extrinsic_params_set_absolute;
    points_relative = points_absolute;

    return (relative_generator_(tracks,
                                camera_views,
                                extrinsic_params_set_relative,
                                points_relative,
                                rotation_similar,
                                translate_similar,
                                scale_similar,
                                camera_id_identity,
                                camera_id_relative));
  }

  Err GenerateNoisedKeysets(const KeysetContainer& keysets_true,
                            KeysetContainer& keysets_noised) const
  {
    typedef EIGEN_VECTOR(Scalar, 2) Key;
    keysets_noised = keysets_true;
    EIGEN_MATRIX(Scalar, 2, 2) covariance_key;
    covariance_key.setIdentity();
    covariance_key *= key_stddev_ * key_stddev_;
    Key min_key;
    min_key << -Scalar(scene_generator_.image_width()) / 2,
               -Scalar(scene_generator_.image_height()) / 2;
    Key max_key = -min_key;
    size_t number_of_keysets = keysets_true.size();
    for (size_t i = 0; i < number_of_keysets; i++)
    {
      size_t number_of_keys = keysets_true[i].size();
      for (size_t j = 0; j < number_of_keys; j++)
      {
        Scalar random;
        hs::math::random::UniformRandomVar<Scalar, 1>::Generate(
          0, 1, random);
        if (random < outlier_ratio_)
        {
          hs::math::random::UniformRandomVar<Scalar, 2>::Generate(
            min_key, max_key, keysets_noised[i][j]);
        }
        else
        {
          Key mean = keysets_true[i][j];
          hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
            mean, covariance_key, keysets_noised[i][j]);
        }
      }
    }

    return 0;
  }

  Err GenerateInitialScene(
    const KeysetContainer& keysets_noised,
    const IntrinsicParamsContainer& intrinsic_params_set,
    const hs::sfm::TrackContainer& tracks,
    const hs::sfm::CameraViewContainer& camera_views,
    const ExtrinsicParamsContainer& extrinsic_params_set_relative,
    const Point3DContainer& points_relative,
    size_t camera_id_identity,
    size_t camera_id_relative,
    ExtrinsicParamsContainer& extrinsic_params_set_relative_estimate,
    Point3DContainer& points_relative_estimate,
    hs::sfm::ObjectIndexMap& image_extrinsic_map,
    hs::sfm::ObjectIndexMap& track_point_map) const
  {
    extrinsic_params_set_relative_estimate.clear();
    extrinsic_params_set_relative_estimate.push_back(
      extrinsic_params_set_relative[camera_id_identity]);
    extrinsic_params_set_relative_estimate.push_back(
      extrinsic_params_set_relative[camera_id_relative]);

    image_extrinsic_map.Resize(extrinsic_params_set_relative.size());
    image_extrinsic_map[camera_id_identity] = 0;
    image_extrinsic_map[camera_id_relative] = 1;

    points_relative_estimate.clear();
    track_point_map.Resize(points_relative.size());
    const hs::sfm::CameraView& view_identity =
      camera_views[camera_id_identity];
    const hs::sfm::CameraView& view_relative =
      camera_views[camera_id_relative];
    for (size_t i = 0; i < view_identity.size(); i++)
    {
      size_t track_id = view_identity[i].first;
      for (size_t j = 0; j < view_relative.size(); j++)
      {
        if (track_id == view_relative[j].first)
        {
          track_point_map[track_id] = points_relative_estimate.size();
          points_relative_estimate.push_back(points_relative[track_id]);
          break;
        }
      }
    }

    return 0;
  }

  Scalar key_stddev() const
  {
    return key_stddev_;
  }

private:
  SceneGenerator scene_generator_;
  KeysetGenerator keyset_generator_;
  RelativeGenerator relative_generator_;
  Scalar outlier_ratio_;
  Scalar key_stddev_;
};

}
}
}

#endif
