#include <iostream>

#include <gtest/gtest.h>

#include "hs_sfm/sfm_utility/synthetic_scene_generator.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"

#include "hs_sfm/incremental/scene_spanner.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestSceneSpanner
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
  typedef hs::sfm::SceneGenerator<Scalar, ImageDimension> SceneGenerator;
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

  typedef hs::sfm::KeysGenerator<Scalar, ImageDimension> KeysGenerator;
  typedef typename KeysGenerator::Keys Keys;
  typedef typename KeysGenerator::KeysContainer KeysContainer;

  typedef hs::sfm::incremental::SceneSpanner<Scalar> Spanner;

public:
  TestSceneSpanner(
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
                     camera_plannar_stddev,
                     camera_rotation_stddev,
                     north_west_angle),
    keys_generator_(image_width, image_height) {}

  Err operator() ()
  {
    //生成相机参数和三维点
    IntrinsicParamsContainer intrinsic_params_set;
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    Point3DContainer points;
    if (scene_generator_(intrinsic_params_set, extrinsic_params_set,
                         images, points) != 0)
    {
      std::cout<<"scene generator failed!\n";
      return -1;
    }

    //生成特征点
    KeysContainer keys_set;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;
    if (keys_generator_(intrinsic_params_set,
                        extrinsic_params_set,
                        points,
                        keys_set,
                        tracks,
                        camera_views) != 0)
    {
      std::cout<<"keys generator failed!\n";
      return -1;
    }

    //计算由相对坐标系转换到绝对坐标系的相似变换七参数
    typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;
    RMatrix rotation_absolute_0 = extrinsic_params_set[0].rotation();
    //rotation_absolute_0.row(0) *= -1;
    //rotation_absolute_0.row(1) *= -1;
    RMatrix rotation = rotation_absolute_0.transpose();
    Position translate = extrinsic_params_set[0].position();
    Scalar scale = (extrinsic_params_set[1].position() - translate).norm();

    //计算真实相对外参数
    RMatrix rotation_absolute_1 = extrinsic_params_set[1].rotation();
    //rotation_absolute_1.row(0) *= -1;
    //rotation_absolute_1.row(1) *= -1;
    const Position& position_absolute = extrinsic_params_set[1].position();
    RMatrix rotation_relative = rotation_absolute_1 * rotation;
    Position position_relative = (position_absolute - translate) / scale;
    position_relative = rotation.transpose() * position_relative;

    //初始化两张影像的外方位元素

    return 0;
  }

private:
  SceneGenerator scene_generator_;
  KeysGenerator keys_generator_;
};

TEST(TestSceneSpanner, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestSceneSpanner<Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 7;
  size_t number_of_cameras_in_strips = 20;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 3000;
  Scalar lateral_overlap_ratio = 0.7;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 100;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 5;
  Scalar north_west_angle = 60;

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
            north_west_angle);

  ASSERT_EQ(0, test());
}

}