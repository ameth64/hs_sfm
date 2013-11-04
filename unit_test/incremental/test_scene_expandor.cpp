#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_sfm/synthetic/keyset_generator.hpp"
#include "hs_sfm/synthetic/relative_generator.hpp"
#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"

#include "hs_sfm/incremental/scene_expandor.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestSceneExpandor
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
  typedef hs::sfm::synthetic::SceneGenerator<Scalar, ImageDimension>
          SceneGenerator;
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

  typedef hs::sfm::synthetic::KeysetGenerator<Scalar, ImageDimension>
          KeysetGenerator;
  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

  typedef hs::sfm::synthetic::RelativeGenerator<Scalar> RelativeGenerator;
  typedef typename RelativeGenerator::RMatrix RMatrix;
  typedef typename RelativeGenerator::Translate Translate;

  typedef hs::sfm::incremental::SceneExpandor<Scalar> Expandor;
  typedef typename Expandor::TrackPointMap TrackPointMap;
  typedef typename Expandor::ImageExtrinsicMap ImageExtrinsicMap;

  typedef hs::sfm::CameraFunctions<Scalar> CameraFunctions;
  typedef typename CameraFunctions::ProjectionMatrix PMatrix;

  typedef EIGEN_VECTOR(Scalar, 2) Key;

public:
  TestSceneExpandor(
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

  Err Test ()
  {
    IntrinsicParamsContainer intrinsic_params_set;
    ExtrinsicParamsContainer extrinsic_params_set_absolute;
    ExtrinsicParamsContainer extrinsic_params_set_relative;
    ImageContainer images;
    Point3DContainer points_absolute;
    Point3DContainer points_relative;
    KeysetContainer keysets_true;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;
    RMatrix rotation_similar;
    Translate translate_similar;
    Scalar scale_similar;
    size_t camera_id_identity;
    size_t camera_id_relative;
    if (TestGenerator(intrinsic_params_set,
                      extrinsic_params_set_absolute,
                      extrinsic_params_set_relative,
                      images,
                      points_absolute,
                      points_relative,
                      keysets_true,
                      tracks,
                      camera_views,
                      rotation_similar,
                      translate_similar,
                      scale_similar,
                      camera_id_identity,
                      camera_id_relative) != 0)
    {
      return -1;
    }

    KeysetContainer keysets_noised;
    if (TestNoiser(keysets_true, keysets_noised) != 0) return -1;

    ExtrinsicParamsContainer extrinsic_params_set_relative_estimate;
    Point3DContainer points_relative_estimate;
    ImageExtrinsicMap image_extrinsic_map;
    TrackPointMap track_point_map;
    if (TestExpandor(keysets_noised,
                    intrinsic_params_set,
                    tracks,
                    camera_views,
                    extrinsic_params_set_relative,
                    points_relative,
                    camera_id_identity,
                    camera_id_relative,
                    extrinsic_params_set_relative_estimate,
                    points_relative_estimate,
                    image_extrinsic_map,
                    track_point_map) != 0)
    {
      return -1;
    }

    if (TestReprojectErr(keysets_noised,
                         intrinsic_params_set,
                         tracks,
                         extrinsic_params_set_relative_estimate,
                         points_relative_estimate,
                         image_extrinsic_map,
                         track_point_map) != 0)
    {
      return -1;
    }

    return 0;
  }

private:
  Err TestGenerator(IntrinsicParamsContainer& intrinsic_params_set,
                    ExtrinsicParamsContainer& extrinsic_params_set_absolute,
                    ExtrinsicParamsContainer& extrinsic_params_set_relative,
                    ImageContainer& images,
                    Point3DContainer& points_absolute,
                    Point3DContainer& points_relative,
                    KeysetContainer& keysets,
                    hs::sfm::TrackContainer& tracks,
                    hs::sfm::CameraViewContainer& camera_views,
                    RMatrix& rotation_similar,
                    Translate& translate_similar,
                    Scalar& scale_similar,
                    size_t& camera_id_identity,
                    size_t& camera_id_relative) const
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

    extrinsic_params_set_relative = extrinsic_params_set_absolute;
    points_relative = points_absolute;
    if (relative_generator_(tracks,
                            camera_views,
                            extrinsic_params_set_relative,
                            points_relative,
                            rotation_similar,
                            translate_similar,
                            scale_similar,
                            camera_id_identity,
                            camera_id_relative) != 0)
    {
      return -1;
    }

    return 0;
  }

  Err TestNoiser(const KeysetContainer& keysets_true,
                 KeysetContainer& keysets_noised) const
  {
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

  Err TestExpandor(
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
    ImageExtrinsicMap& image_extrinsic_map,
    TrackPointMap& track_point_map) const
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
    const hs::sfm::CameraView& view_identity = camera_views[camera_id_identity];
    const hs::sfm::CameraView& view_relative = camera_views[camera_id_relative];
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

    Expandor expandor(50, key_stddev_ * 4);
    if (expandor(keysets_noised, intrinsic_params_set, tracks,
                 extrinsic_params_set_relative_estimate,
                 image_extrinsic_map,
                 points_relative_estimate,
                 track_point_map) != 0)
    {
      return -1;
    }

    return 0;
  }

  Err TestReprojectErr(
    const KeysetContainer& keysets_noised,
    const IntrinsicParamsContainer& intrinsic_params_set,
    const hs::sfm::TrackContainer& tracks,
    const ExtrinsicParamsContainer& extrinsic_params_set_relative_estimate,
    const Point3DContainer& points_relative_estimate,
    const ImageExtrinsicMap& image_extrinsic_map,
    const TrackPointMap& track_point_map) const
  {
    auto itr_track_point = track_point_map.begin();
    auto itr_track_point_end = track_point_map.end();
    size_t number_of_reprojections = 0;
    Scalar mean_reprojection_error = 0;
    for (; itr_track_point != itr_track_point_end; ++itr_track_point)
    {
      size_t track_id = itr_track_point->first;
      size_t point_id = itr_track_point->second;
      const Point3D& point = points_relative_estimate[point_id];
      const hs::sfm::Track& track  = tracks[track_id];
      size_t number_of_views = track.size();
      for (size_t i = 0; i < number_of_views; i++)
      {
        size_t image_id = track[i].first;
        size_t key_id = track[i].second;
        auto itr_image_extrinsic = image_extrinsic_map.find(image_id);
        if (itr_image_extrinsic != image_extrinsic_map.end())
        {
          size_t extrinsic_id = itr_image_extrinsic->second;
          PMatrix P = CameraFunctions::GetProjectionMatrix(
                        intrinsic_params_set[image_id],
                        extrinsic_params_set_relative_estimate[extrinsic_id]);
          Point3D reprojection = P.block(0, 0, 3, 3) * point + P.col(3);
          reprojection /= reprojection[2];
          Scalar reprojection_error = (keysets_noised[image_id][key_id] -
                                       reprojection.segment(0, 2)).norm();
          mean_reprojection_error += reprojection_error;
          number_of_reprojections++;
        }
      }
    }
    mean_reprojection_error /= Scalar(number_of_reprojections);
    std::cout<<"mean reprojection error:"<<mean_reprojection_error<<"\n;";
    if (mean_reprojection_error < 2 * key_stddev_)
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }

private:
  SceneGenerator scene_generator_;
  KeysetGenerator keyset_generator_;
  RelativeGenerator relative_generator_;
  Scalar outlier_ratio_;
  Scalar key_stddev_;
};

TEST(TestSceneExpandor, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestSceneExpandor<Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 5;
  size_t number_of_cameras_in_strips = 8;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 40000;
  Scalar lateral_overlap_ratio = 0.7;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 100;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 5;
  Scalar north_west_angle = 60;
  Scalar outlier_ratio = 0.4;
  Scalar key_stddev = 2.0;

  Test test(focal_length_in_metre,
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

  ASSERT_EQ(0, test.Test());
}

}
