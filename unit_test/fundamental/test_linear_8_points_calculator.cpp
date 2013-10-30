#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/synthetic/keyset_generator.hpp"

#include "hs_sfm/fundamental/linear_8_points_calculator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestLinear8PointsCalculator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
  typedef hs::sfm::synthetic::SceneGenerator<Scalar, ImageDimension>
          SceneGenerator;
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

  typedef hs::sfm::synthetic::KeysetGenerator<Scalar, ImageDimension>
          KeysetGenerator;
  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

  typedef hs::sfm::fundamental::Linear8PointsCalculator<Scalar> Calculator;
  typedef typename Calculator::Key Key;
  typedef typename Calculator::KeyContainer KeyContainer;
  typedef typename Calculator::FMatrix FMatrix;
  typedef EIGEN_VECTOR(Scalar, 3) HKey;
  typedef EIGEN_VECTOR(Scalar, 3) HLine;

public:
  TestLinear8PointsCalculator(
    Scalar focal_length_in_metre,
    size_t number_of_strps,
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
                       number_of_strps > 1 ? 1 : 2,
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
      keys_generator_(image_width, image_height){}

  Err operator()()
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
    KeysetContainer keysets;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;
    if (keys_generator_(intrinsic_params_set,
                        extrinsic_params_set,
                        points,
                        keysets,
                        tracks,
                        camera_views) != 0)
    {
      std::cout<<"keys generator failed!\n";
      return -1;
    }

    //生成计算F矩阵所需的匹配点
    KeyContainer keys_left;
    KeyContainer keys_right;
    size_t number_of_tracks = tracks.size();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      if (tracks[i].size() == 2)
      {
        size_t key_left_id = tracks[i][0].second;
        size_t key_right_id = tracks[i][1].second;
        Key key_left = keysets[0][key_left_id];
        Key key_right = keysets[1][key_right_id];
        keys_left.push_back(key_left);
        keys_right.push_back(key_right);
      }
    }

    //为匹配点加入高斯误差
    size_t number_of_keys = keys_left.size();
    EIGEN_MATRIX(Scalar, 2, 2) key_covariance;
    key_covariance.setIdentity();
    Scalar key_stddev = 2;
    key_covariance *= key_stddev * key_stddev;
    for (size_t i = 0; i < number_of_keys; i++)
    {
      Key key_left_mean = keys_left[i];
      hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
        key_left_mean, key_covariance, keys_left[i]);
      Key key_right_mean = keys_right[i];
      hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
        key_right_mean, key_covariance, keys_right[i]);
    }

    //计算F矩阵
    Calculator calculator;
    FMatrix f_matrix;
    if (calculator(keys_left, keys_right, f_matrix) != 0)
    {
      std::cout<<"calculator failed!\n";
      return -1;
    }

    //检查对极线误差
    Scalar mean_error = Scalar(0);
    for (size_t i = 0; i < number_of_keys; i++)
    {
      const Key& key_left = keys_left[i];
      const Key& key_right = keys_right[i];
      HKey hkey_left;
      hkey_left << key_left[0],
                   key_left[1],
                   1;
      HKey hkey_right;
      hkey_right << key_right[0],
                    key_right[1],
                    1;
      HLine epiline_left = f_matrix * hkey_left;
      HLine epiline_right = f_matrix.transpose() * hkey_right;
      Scalar error_left = std::abs(hkey_right.dot(epiline_left)) /
                          epiline_left.segment(0, 2).norm();
      Scalar error_right = std::abs(hkey_left.dot(epiline_right)) /
                           epiline_right.segment(0, 2).norm();
      mean_error += error_left + error_right;
    }
    mean_error /= Scalar(number_of_keys);

    //TODO:该阈值设置较随意！应改善！
    if (mean_error <= key_stddev * Scalar(2) + 1)
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
  KeysetGenerator keys_generator_;

};

TEST(TestLinear8PointsCalculator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestLinear8PointsCalculator<Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 1;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 10000;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_plannar_stddev = 5;
  Scalar camera_rot_stddev = 1;
  Scalar north_west_angle = 60;

  Test test(focal_length_in_metre,
            number_of_strips,
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
            north_west_angle);

  ASSERT_EQ(0, test());
}

}