#include <iostream>
#include <set>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_sfm/synthetic/flight_generator.hpp"
#include "hs_sfm/synthetic/multiple_camera_keyset_generator.hpp"
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
  typedef hs::sfm::synthetic::FlightGenerator<Scalar, ImageDimension>
          FlightGenerator;
  typedef typename FlightGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename FlightGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename FlightGenerator::Image Image;
  typedef typename FlightGenerator::ImageContainer ImageContainer;
  typedef typename FlightGenerator::Point3D Point3D;
  typedef typename FlightGenerator::Point3DContainer Point3DContainer;

  typedef hs::sfm::synthetic::MultipleCameraKeysetGenerator<Scalar,
                                                            ImageDimension>
          KeysetGenerator;
  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

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
  typedef typename KeysetGenerator::IntrinsicParams IntrinsicParams;
  typedef typename KeysetGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;

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
    Scalar outlier_ratio,
    const IntrinsicParamsContainer& intrinsic_params_set)
    : flight_generator_(focal_length_in_metre,
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
                       camera_rot_stddev),
      outlier_ratio_(outlier_ratio),
      intrinsic_params_set_(intrinsic_params_set),
      image_width_(image_width),
      image_height_(image_height){}

  Err operator()()
  {
    //生成相机参数和三维点
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    Point3DContainer points;
    if (flight_generator_(extrinsic_params_set, images, points) != 0)
    {
      std::cout<<"scene generator failed!\n";
      return -1;
    }

    //生成特征点
    std::vector<size_t> image_intrinsic_map;
    image_intrinsic_map.push_back(0);
    image_intrinsic_map.push_back(1);
    KeysetContainer keysets;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;
    if (keys_generator_(intrinsic_params_set_,
                        extrinsic_params_set,
                        images,
                        points,
                        image_intrinsic_map,
                        keysets,
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
        key_pair.first = keysets[0][key_left_id];
        key_pair.second = keysets[1][key_right_id];
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
    max << Scalar(flight_generator_.image_width()),
           Scalar(flight_generator_.image_height());
    Key min;
    min.setZero();
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
    if (refiner(key_pairs, key_stddev * 2 * 6,
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
  FlightGenerator flight_generator_;
  KeysetGenerator keys_generator_;
  IntrinsicParamsContainer intrinsic_params_set_;
  Scalar outlier_ratio_;
  ImageDimension image_width_;
  ImageDimension image_height_;
};

TEST(TestLinear8PointsRansacRefiner, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestLinear8PointsRansacRefiner<Scalar, ImageDimension> Test;
  typedef Test::IntrinsicParams IntrinsicParams;
  typedef Test::IntrinsicParamsContainer IntrinsicParamsContainer;

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

  IntrinsicParams intrinsic_params_0(4835.47665904517026,
                                     0,
                                     3000-42.4095312016,
                                     2000+31.699212823,
                                     1,
                                     -0.0050490462006048831,
                                     0.031293804298609881,
                                     -0.030794820960442223,
                                     -0.00055376548320189127,
                                     -0.00049877717768381476);
  IntrinsicParams intrinsic_params_1(4886.17891666633641,
                                     0,
                                     3000-35.2052431556,
                                     2000+16.4262220759,
                                     1,
                                     -0.10316088386868619,
                                     0.13520490482776426,
                                     -0.05489235547426094,
                                     4.1434720317373253e-006,
                                     -0.00025018439997095336);
  IntrinsicParamsContainer intrinsic_params_set;
  intrinsic_params_set.push_back(intrinsic_params_0);
  intrinsic_params_set.push_back(intrinsic_params_1);

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
            outlier_ration,
            intrinsic_params_set);

  ASSERT_EQ(0, test());
}

}
