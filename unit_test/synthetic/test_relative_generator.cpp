#include <iostream>
#include <cmath>

#include <gtest/gtest.h>

#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/synthetic/keyset_generator.hpp"

#include "hs_sfm/synthetic/relative_generator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestRelativeGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;
  typedef hs::sfm::synthetic::SceneGenerator<Scalar, ImageDimension>
          SceneGenerator;
  typedef hs::sfm::synthetic::KeysetGenerator<Scalar, ImageDimension>
          KeysetGenerator;
  typedef hs::sfm::synthetic::RelativeGenerator<Scalar>
          RelativeGenerator;

private:
  typedef typename SceneGenerator::IntrinsicParams IntrinsicParams;
  typedef typename IntrinsicParams::KMatrix KMatrix;
  typedef typename SceneGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename ExtrinsicParams::Position Position;
  typedef typename SceneGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename SceneGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename SceneGenerator::Point3D Point3D;
  typedef typename SceneGenerator::Point3DContainer Point3DContainer;
  typedef typename SceneGenerator::Image Image;
  typedef typename SceneGenerator::ImageContainer ImageContainer;

  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

  typedef typename RelativeGenerator::RMatrix RMatrix;
  typedef typename RelativeGenerator::Translate Translate;

public:
  TestRelativeGenerator(
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
    Scalar camera_planar_stddev,
    Scalar camera_rotation_stddev,
    Scalar north_west_angle)
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
                       camera_planar_stddev,
                       camera_rotation_stddev,
                       north_west_angle),
      keyset_generator_(image_width, image_height) {}

  Err Test() const
  {
    IntrinsicParamsContainer intrinsic_params_set;
    ExtrinsicParamsContainer extrinsic_params_set_absolute;
    ExtrinsicParamsContainer extrinsic_params_set_relative;
    ImageContainer images;
    Point3DContainer points_absolute;
    Point3DContainer points_relative;
    KeysetContainer keysets;
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
                      keysets,
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

    if (TestPair(extrinsic_params_set_relative,
                 camera_id_identity,
                 camera_id_relative) != 0)
    {
      return -1;
    }

    if (TestReprojectError(intrinsic_params_set,
                           extrinsic_params_set_relative,
                           images,
                           points_relative,
                           keysets,
                           tracks,
                           camera_views) != 0)
    {
      return -1;
    }

    if (TestSimilar(extrinsic_params_set_relative,
                    extrinsic_params_set_absolute,
                    points_relative,
                    points_absolute,
                    rotation_similar,
                    translate_similar,
                    scale_similar) != 0)
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

  Err TestPair(const ExtrinsicParamsContainer& extrinsic_params_set_relative,
               Scalar camera_id_identity,
               Scalar camera_id_relative) const
  {
    Scalar threshold = 1e-8;
    RMatrix rotation_identity =
      extrinsic_params_set_relative[camera_id_identity].rotation();
    Position position_identity =
      extrinsic_params_set_relative[camera_id_identity].position();

    if (!rotation_identity.isApprox(RMatrix::Identity(), threshold))
    {
      return -1;
    }

    if (!position_identity.isApprox(Position::Zero(), threshold))
    {
      return -1;
    }

    Position position_relative =
      extrinsic_params_set_relative[camera_id_relative].position();
    if (std::abs(position_relative.norm() - Scalar(1)) > threshold)
    {
      return -1;
    }

    return 0;
  }

  Err TestReprojectError(
    const IntrinsicParamsContainer& intrinsic_params_set,
    const ExtrinsicParamsContainer& extrinsic_params_set_relative,
    const ImageContainer& images,
    const Point3DContainer& points_relative,
    const KeysetContainer& keysets,
    const hs::sfm::TrackContainer& tracks,
    const hs::sfm::CameraViewContainer& camera_views) const
  {
    typedef EIGEN_MATRIX(Scalar, 3, 4) PMatrix;
    size_t number_of_tracks = tracks.size();
    size_t number_of_track_reprojections = 0;
    Scalar threshold = 1e-8;
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();
      const Point3D& point = points_relative[i];
      for (size_t j = 0; j < number_of_views; j++)
      {
        size_t camera_id = tracks[i][j].first;
        size_t key_id = tracks[i][j].second;
        RMatrix R = extrinsic_params_set_relative[camera_id].rotation();
        Position c = extrinsic_params_set_relative[camera_id].position();
        Position t = -R * c;
        KMatrix K = intrinsic_params_set[camera_id].GetKMatrix();
        PMatrix P;
        P.block(0, 0, 3, 3) = R;
        P.block(0, 3, 3, 1) = t;
        P = K * P;
        Point3D hkey = P.block(0, 0, 3, 3) * point + P.col(3);
        hkey /= hkey[2];
        EIGEN_VECTOR(Scalar, 2) key = hkey.segment(0, 2);

        if (!key.isApprox(keysets[camera_id][key_id], threshold))
        {
          return -1;
        }

        number_of_track_reprojections++;
      }
    }

    size_t number_of_cameras = camera_views.size();
    size_t number_of_camera_views_reprojections = 0;
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      RMatrix R = extrinsic_params_set_relative[i].rotation();
      Position c = extrinsic_params_set_relative[i].position();
      Position t = -R * c;
      KMatrix K = intrinsic_params_set[i].GetKMatrix();
      PMatrix P;
      P.block(0, 0, 3, 3) = R;
      P.block(0, 3, 3, 1) = t;
      P = K * P;
      size_t number_of_views = camera_views[i].size();
      for (size_t j = 0; j < number_of_views; j++)
      {
        size_t point_id = camera_views[i][j].first;
        size_t key_id = camera_views[i][j].second;
        const Point3D point = points_relative[point_id];
        Point3D hkey = P.block(0, 0, 3, 3) * point + P.col(3);
        hkey /= hkey[2];
        EIGEN_VECTOR(Scalar, 2) key = hkey.segment(0, 2);
        if (!key.isApprox(keysets[i][key_id], threshold))
        {
          return -1;
        }

        number_of_camera_views_reprojections++;
      }
    }

    if (number_of_track_reprojections != number_of_camera_views_reprojections)
    {
      return -1;
    }

    return 0;
  }

  Err TestSimilar(const ExtrinsicParamsContainer& extrinsic_params_set_relative,
                  const ExtrinsicParamsContainer& extrinsic_params_set_absolute,
                  const Point3DContainer& points_relative,
                  const Point3DContainer& points_absolute,
                  const RMatrix& rotation_similar,
                  const Translate& translate_similar,
                  Scalar scale_similar) const
  {
    Scalar threshold = 1e-8;
    size_t number_of_cameras = extrinsic_params_set_absolute.size();
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      RMatrix rotation_relative =
        extrinsic_params_set_relative[i].rotation();
      Position position_relative =
        extrinsic_params_set_relative[i].position();
      RMatrix rotation_absolute =
        extrinsic_params_set_absolute[i].rotation();
      Position position_absolute =
        extrinsic_params_set_absolute[i].position();

      RMatrix rotation_absolute_estimate =
        rotation_relative * rotation_similar.transpose();
      Position position_absolute_estimate =
        scale_similar * rotation_similar * position_relative +
        translate_similar;

      if (!rotation_absolute_estimate.isApprox(rotation_absolute, threshold))
      {
        return -1;
      }

      if (!position_absolute_estimate.isApprox(position_absolute, threshold))
      {
        return -1;
      }
    }

    size_t number_of_points = points_absolute.size();
    for (size_t i = 0; i < number_of_points; i++)
    {
      Point3D point_relative = points_relative[i];
      Point3D point_absolute = points_absolute[i];
      Point3D point_absolute_estimate =
        scale_similar * rotation_similar * point_relative + translate_similar;
      if (!point_absolute_estimate.isApprox(point_absolute, threshold))
      {
        return -1;
      }
    }

    return 0;
  }

private:
  SceneGenerator scene_generator_;
  KeysetGenerator keyset_generator_;
  RelativeGenerator relative_generator_;
};

TEST(TestRelativeGenerator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestRelativeGenerator<Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 15;
  size_t number_of_cameras_in_strip = 20;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 2000;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 2;
  Scalar camera_planar_stddev = 2;
  Scalar camera_rotation_stddev = 10;
  Scalar north_west_angle = 60;

  Test test(focal_length_in_metre,
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
            camera_planar_stddev,
            camera_rotation_stddev,
            north_west_angle);

  ASSERT_EQ(0, test.Test());
}

}
