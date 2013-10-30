#include <iostream>
#include <set>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/synthetic/keyset_generator.hpp"

#include "hs_sfm/essential/ematrix_5_points_ransac_refiner.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestEMatrix5PointsRansacRefiner
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

  typedef hs::sfm::essential::EMatrix5PointsRansacRefiner<Scalar> Refiner;
  typedef typename Refiner::HKey HKey;
  typedef typename Refiner::HKeyPair HKeyPair;
  typedef typename Refiner::HKeyPairContainer HKeyPairContainer;
  typedef typename Refiner::IndexSet IndexSet;

public:
  TestEMatrix5PointsRansacRefiner(
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

    //生成计算E矩阵所需的匹配点，并加入高斯误差和野值
    HKeyPairContainer key_pairs;
    size_t number_of_tracks = tracks.size();
    EIGEN_MATRIX(Scalar, 2, 2) key_covariance;
    key_covariance.setIdentity();
    Scalar key_stddev = 2;
    key_covariance *= key_stddev * key_stddev;
    KMatrix K_left_inverse = intrinsic_params_set[0].GetKMatrix().inverse();
    KMatrix K_right_inverse = intrinsic_params_set[1].GetKMatrix().inverse();
    EIGEN_VECTOR(Scalar, 2) max;
    max << Scalar(image_width_) / (Scalar(2)),
           Scalar(image_height_) / (Scalar(2));
    EIGEN_VECTOR(Scalar, 2) min = -max;
    std::set<size_t> true_outlier_indices;
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      if (tracks[i].size() == 2)
      {
        EIGEN_VECTOR(Scalar, 2) noised_key_left;
        EIGEN_VECTOR(Scalar, 2) noised_key_right;
        size_t key_left_id = tracks[i][0].second;
        size_t key_right_id = tracks[i][1].second;
        HKeyPair key_pair;
        hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
          keysets[0][key_left_id], key_covariance, noised_key_left);
        hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
          keysets[1][key_right_id], key_covariance, noised_key_right);

        Scalar random;
        hs::math::random::UniformRandomVar<Scalar, 1>::Generate(
          Scalar(0), Scalar(1), random);
        if (random < outlier_ratio_)
        {
          hs::math::random::UniformRandomVar<Scalar, 2>::Generate(
            min, max, noised_key_left);
          hs::math::random::UniformRandomVar<Scalar, 2>::Generate(
            min, max, noised_key_right);
          true_outlier_indices.insert(key_pairs.size());
        }

        key_pair.first.template segment<2>(0) = noised_key_left;
        key_pair.second.template segment<2>(0) = noised_key_right;
        key_pair.first[2] = Scalar(1);
        key_pair.second[2] = Scalar(1);
        key_pair.first = K_left_inverse * key_pair.first;
        key_pair.second = K_right_inverse * key_pair.second;
        key_pairs.push_back(key_pair);
      }
    }

    Refiner refiner;
    HKeyPairContainer refined_key_pairs;
    IndexSet estimated_inlier_indices;
    if (refiner(key_pairs,
                key_stddev * 2 * 4 / scene_generator_.GetFocalLengthInPixel(),
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
  KeysetGenerator keys_generator_;
  Scalar outlier_ratio_;
  ImageDimension image_width_;
  ImageDimension image_height_;
};

TEST(TestEMatrix5PointsRansacRefiner, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestEMatrix5PointsRansacRefiner<Scalar, ImageDimension> Test;

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