#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_sfm/sfm_utility/synthetic_scene_generator.hpp"

#include "hs_sfm/projective/pmatrix_dlt_ransac_refiner.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestPMatrixDLTRansacRefiner
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

  typedef hs::sfm::projective::PMatrixDLTRansacRefiner<Scalar> Refiner;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef typename Refiner::Correspondence Correspondence;
  typedef typename Refiner::CorrespondenceContainer CorrespondenceContainer;
  typedef typename Refiner::IndexSet IndexSet;

public:
  TestPMatrixDLTRansacRefiner(
    Scalar focal_length_in_metre,
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    size_t number_of_points,
    Scalar scene_max_height,
    Scalar camera_height_stddev,
    Scalar camera_planar_stddev,
    Scalar camera_rot_stdddev,
    Scalar north_west_angle,
    Scalar outlier_ratio)
    : scene_generator_(focal_length_in_metre,
                       1,
                       1,
                       ground_resolution,
                       image_width,
                       image_height,
                       pixel_size,
                       number_of_points,
                       0.6,
                       0.8,
                       scene_max_height,
                       camera_height_stddev,
                       camera_planar_stddev,
                       camera_rot_stdddev,
                       north_west_angle),
      keys_generator_(image_width, image_height),
      outlier_ratio_(outlier_ratio),
      image_width_(image_width),
      image_height_(image_height) {}

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

    //生成计算P矩阵所需的对应点
    CorrespondenceContainer correspondences;
    size_t number_of_tracks = tracks.size();
    Scalar key_stddev = 2;
    EIGEN_MATRIX(Scalar, 2, 2) key_covariance;
    key_covariance.setIdentity();
    key_covariance *= key_stddev * key_stddev;
    EIGEN_VECTOR(Scalar, 2) max;
    max << Scalar(image_width_) / (Scalar(2)),
           Scalar(image_height_) / (Scalar(2));
    EIGEN_VECTOR(Scalar, 2) min = -max;
    std::set<size_t> true_outlier_indices;
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      if (tracks[i].size() == 1)
      {
        Key key_mean = keys_set[0][tracks[i][0].second];
        Point point_mean = points[i];
        Correspondence correspondence(key_mean, point_mean);
        hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
          key_mean, key_covariance, correspondence.first);
        Scalar random;
        hs::math::random::UniformRandomVar<Scalar, 1>::Generate(
          Scalar(0), Scalar(1), random);
        if (random < outlier_ratio_)
        {
          hs::math::random::UniformRandomVar<Scalar, 2>::Generate(
            min, max, correspondence.first);
          true_outlier_indices.insert(correspondences.size());
        }
        correspondences.push_back(correspondence);
      }
    }

    Refiner refiner;
    CorrespondenceContainer refined_correspondences;
    IndexSet estimated_inlier_indices;
    if (refiner(correspondences,
                key_stddev * 4,
                refined_correspondences,
                estimated_inlier_indices) != 0)
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

TEST(TestPMatrixDLTRansacRefiner, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestPMatrixDLTRansacRefiner<Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 5000;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_plannar_stddev = 5;
  Scalar camera_rot_stddev = 1;
  Scalar north_west_angle = 60;
  Scalar outlier_ration = 0.2;

  Test test(focal_length_in_metre,
            ground_resolution,
            image_width,
            image_height,
            pixel_size,
            number_of_points,
            scene_max_height,
            camera_height_stddev,
            camera_plannar_stddev,
            camera_rot_stddev,
            north_west_angle,
            outlier_ration);

  ASSERT_EQ(0, test());

}

}