#include <iostream>
#include <set>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_sfm/sfm_utility/synthetic_scene_generator.hpp"
#include "hs_sfm/fundamental/linear_8_points_calculator.hpp"

#include "hs_sfm/fundamental/linear_8_points_ransac_refiner.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestLinear8PointsRansacRefiner
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
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

  typedef hs::sfm::fundamental::Linear8PointsCalculator<Scalar> Calculator;
  typedef typename Calculator::Key Key;
  typedef typename Calculator::FMatrix FMatrix;
  typedef EIGEN_VECTOR(Scalar, 3) HKey;
  typedef EIGEN_VECTOR(Scalar, 3) HLine;

  typedef hs::sfm::fundamental::Linear8PointsRansacRefiner<Scalar> Refiner;
  typedef typename Refiner::KeyPair KeyPair;
  typedef typename Refiner::KeyPairContainer KeyPairContainer;
  typedef typename Refiner::IndexSet IndexSet;

public:
  TestLinear8PointsRansacRefiner(
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
    Scalar north_west_angle,
    Scalar outlier_ratio)
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
      keys_generator_(image_width, image_height),
      outlier_ratio_(outlier_ratio),
      image_width_(image_width),
      image_height_(image_height){}

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

    //生成计算F矩阵所需的匹配点
    KeyPairContainer key_pairs;
    size_t number_of_tracks = tracks.size();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      if (tracks[i].size() == 2)
      {
        size_t key_left_id = tracks[i][0].second;
        size_t key_right_id = tracks[i][1].second;
        KeyPair key_pair;
        key_pair.first = keys_set[0][key_left_id];
        key_pair.second = keys_set[1][key_right_id];
        key_pairs.push_back(key_pair);
      }
    }

    //为匹配点加入高斯误差
    size_t number_of_keys = key_pairs.size();
    EIGEN_MATRIX(Scalar, 2, 2) key_covariance;
    key_covariance.setIdentity();
    Scalar key_stddev = 2;
    key_covariance *= key_stddev * key_stddev;
    for (size_t i = 0; i < number_of_keys; i++)
    {
      Key key_left_mean = key_pairs[i].first;
      hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
        key_left_mean, key_covariance, key_pairs[i].first);
      Key key_right_mean = key_pairs[i].second;
      hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
        key_right_mean, key_covariance, key_pairs[i].second);
    }

    //为匹配点加入野值
    Key max;
    max << Scalar(image_width_) / (Scalar(2)),
           Scalar(image_height_) / (Scalar(2));
    Key min = -max;
    std::set<size_t> true_outlier_indices;
    for (size_t i = 0; i < number_of_keys; i++)
    {
      Scalar random;
      hs::math::random::UniformRandomVar<Scalar, 1>::Generate(
        Scalar(0), Scalar(1), random);
      if (random < outlier_ratio_)
      {
        hs::math::random::UniformRandomVar<Scalar, 2>::Generate(
          min, max, key_pairs[i].first);
        hs::math::random::UniformRandomVar<Scalar, 2>::Generate(
          min, max, key_pairs[i].second);
        true_outlier_indices.insert(i);
      }
    }

    Refiner refiner;
    KeyPairContainer refined_key_pairs;
    IndexSet estimated_inlier_indices;
    if (refiner(key_pairs, key_stddev * 2,
                refined_key_pairs, estimated_inlier_indices) != 0)
    {
      std::cout<<"refiner failed!\n";
      return -1;
    }

    //统计被错误作为内点的数量
    size_t number_of_estimated_inliers = estimated_inlier_indices.size();
    size_t number_of_false_inliers = 0;
    for (size_t i = 0; i < number_of_estimated_inliers; i++)
    {
      if (true_outlier_indices.find(estimated_inlier_indices[i]) !=
          true_outlier_indices.end())
      {
        number_of_false_inliers++;
      }
    }

    if (Scalar(number_of_false_inliers) /
        Scalar(number_of_estimated_inliers) < 0.01)
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }

  void set_outlier_ratio(Scalar outlier_ratio)
  {
    outlier_ratio_ = outlier_ratio;
  }

private:
  SceneGenerator scene_generator_;
  KeysGenerator keys_generator_;
  Scalar outlier_ratio_;
  ImageDimension image_width_;
  ImageDimension image_height_;
};

TEST(TestLinear8PointsRansacRefiner, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestLinear8PointsRansacRefiner<Scalar, ImageDimension> Test;

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
  Scalar outlier_ration = 0.2;

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
            north_west_angle,
            outlier_ration);

  ASSERT_EQ(0, test());
}

}