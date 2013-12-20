#include <iostream>

#include <gtest/gtest.h>

#include "hs_sfm/synthetic/scene_generator.hpp"

#include "hs_sfm/synthetic/keyset_generator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestKeysetGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;
  typedef hs::sfm::synthetic::SceneGenerator<Scalar, ImageDimension>
          SceneGenerator;
  typedef hs::sfm::synthetic::KeysetGenerator<Scalar, ImageDimension>
          KeysetGenerator;

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

  typedef hs::sfm::CameraFunctions<Scalar> CameraFunctions;
  typedef typename CameraFunctions::ProjectionMatrix PMatrix;

  typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;

public:
  TestKeysetGenerator(
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
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    Point3DContainer points;
    KeysetContainer keysets;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;
    if (TestGenerator(intrinsic_params_set,
                      extrinsic_params_set,
                      images,
                      points,
                      keysets,
                      tracks,
                      camera_views) != 0)
    {
      return -1;
    }

    if (TestReprojectError(intrinsic_params_set,
                           extrinsic_params_set,
                           images,
                           points,
                           keysets,
                           tracks,
                           camera_views) != 0)
    {
      return -1;
    }

    return 0;
  }

private:
  Err TestGenerator(IntrinsicParamsContainer& intrinsic_params_set,
                    ExtrinsicParamsContainer& extrinsic_params_set,
                    ImageContainer& images,
                    Point3DContainer& points,
                    KeysetContainer& keysets,
                    hs::sfm::TrackContainer& tracks,
                    hs::sfm::CameraViewContainer& camera_views) const
  {
    if (scene_generator_(intrinsic_params_set,
                         extrinsic_params_set,
                         images,
                         points) != 0)
    {
      return -1;
    }

    if (keyset_generator_(intrinsic_params_set,
                          extrinsic_params_set,
                          points,
                          keysets,
                          tracks,
                          camera_views) != 0)
    {
      return -1;
    }

    return 0;
  }

  Err TestReprojectError(const IntrinsicParamsContainer& intrinsic_params_set,
                         const ExtrinsicParamsContainer& extrinsic_params_set,
                         const ImageContainer& images,
                         const Point3DContainer& points,
                         const KeysetContainer& keysets,
                         const hs::sfm::TrackContainer& tracks,
                         const hs::sfm::CameraViewContainer& camera_views) const
  {
    size_t number_of_tracks = tracks.size();
    size_t number_of_track_reprojections = 0;
    Scalar threshold = 1e-8;
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();
      const Point3D& point = points[i];
      for (size_t j = 0; j < number_of_views; j++)
      {
        size_t camera_id = tracks[i][j].first;
        size_t key_id = tracks[i][j].second;
        RMatrix R = extrinsic_params_set[camera_id].rotation();
        Position c = extrinsic_params_set[camera_id].position();
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
      RMatrix R = extrinsic_params_set[i].rotation();
      Position c = extrinsic_params_set[i].position();
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
        const Point3D point = points[point_id];
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

private:
  SceneGenerator scene_generator_;
  KeysetGenerator keyset_generator_;
};

TEST(TestKeysetGenerator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestKeysetGenerator<Scalar, ImageDimension> Test;

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
