#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/synthetic/flight_generator.hpp"
#include "hs_sfm/synthetic/multiple_camera_keyset_generator.hpp"

#include "hs_sfm/essential/ematrix_5_points_calculator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestEMatrix5PointsCalculator
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

  typedef hs::sfm::essential::EMatrix5PointsCalculator<Scalar> Calculator;
  typedef typename Calculator::EMatrix EMatrix;
  typedef typename Calculator::EMatrixHypotheses EMatrixHypotheses;
  typedef typename Calculator::HKey HKey;
  typedef typename Calculator::HKeyPair HKeyPair;
  typedef typename Calculator::HKeyPairContainer HKeyPairContainer;
public:
  typedef typename KeysetGenerator::IntrinsicParams IntrinsicParams;
  typedef typename KeysetGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;

private:
  typedef typename IntrinsicParams::KMatrix KMatrix;
  typedef EIGEN_VECTOR(Scalar, 3) HLine;

public:
  TestEMatrix5PointsCalculator(
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
      intrinsic_params_set_(intrinsic_params_set) {}

  Err operator() ()
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

    //生成计算E矩阵所需的匹配点，并加入高斯误差
    HKeyPairContainer key_pairs;
    size_t number_of_tracks = tracks.size();
    EIGEN_MATRIX(Scalar, 2, 2) key_covariance;
    key_covariance.setIdentity();
    Scalar key_stddev = 2;
    key_covariance *= key_stddev * key_stddev;
    KMatrix K_left_inverse = intrinsic_params_set_[0].GetKMatrix().inverse();
    KMatrix K_right_inverse = intrinsic_params_set_[1].GetKMatrix().inverse();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      if (tracks[i].size() == 2)
      {
        EIGEN_VECTOR(Scalar, 2) noised_key;
        size_t key_left_id = tracks[i][0].second;
        size_t key_right_id = tracks[i][1].second;
        HKeyPair key_pair;
        hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
          keysets[0][key_left_id], key_covariance, noised_key);
          key_pair.first.template segment<2>(0) = noised_key;
        hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
          keysets[1][key_right_id], key_covariance, noised_key);
          key_pair.second.template segment<2>(0) = noised_key;
        key_pair.first[2] = Scalar(1);
        key_pair.second[2] = Scalar(1);
        key_pair.first = K_left_inverse * key_pair.first;
        key_pair.second = K_right_inverse * key_pair.second;
        key_pairs.push_back(key_pair);
      }
    }

    //计算E矩阵
    Calculator calculator;
    EMatrix e_matrix;
    if (calculator(key_pairs, e_matrix) != 0)
    {
      std::cout<<"calculator failed!\n";
      return -1;
    }

    //检查对极线误差
    Scalar mean_error = Scalar(0);
    size_t number_of_key_pairs = key_pairs.size();
    Scalar focal_length_in_pix = flight_generator_.GetFocalLengthInPixel();
    for (size_t i = 0; i < number_of_key_pairs; i++)
    {
      const HKey& hkey_left = key_pairs[i].first;
      const HKey& hkey_right = key_pairs[i].second;

      HLine epiline_left = e_matrix * hkey_left;
      HLine epiline_right = e_matrix.transpose() * hkey_right;
      Scalar error_left = std::abs(hkey_right.dot(epiline_left)) /
        epiline_left.segment(0, 2).norm() *
        intrinsic_params_set_[0].focal_length();
      Scalar error_right = std::abs(hkey_left.dot(epiline_right)) /
        epiline_right.segment(0, 2).norm() *
        intrinsic_params_set_[1].focal_length();
      mean_error += error_left + error_right;
    }
    mean_error /= Scalar(number_of_key_pairs);

    //TODO:该阈值设置较随意！应改善！
    std::cout<<"mean_error:"<<mean_error<<"\n";
    if (mean_error <= 16)
    {
      return 0;
    }
    else
    {
      return -1;
    }

    return 0;
  }
private:
  FlightGenerator flight_generator_;
  KeysetGenerator keys_generator_;
  IntrinsicParamsContainer intrinsic_params_set_;
};

TEST(TestEMatrix5PointsCalculator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestEMatrix5PointsCalculator<Scalar, ImageDimension> Test;
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

  IntrinsicParams intrinsic_params_0(4835.47665904517026,
                                     0,
                                     3000-42.4095312016,
                                     2000+31.699212823,
                                     1/*,
                                     -0.0050490462006048831,
                                     0.031293804298609881,
                                     -0.030794820960442223,
                                     -0.00055376548320189127,
                                     -0.00049877717768381476*/);
  IntrinsicParams intrinsic_params_1(4886.17891666633641,
                                     0,
                                     3000-35.2052431556,
                                     2000+16.4262220759,
                                     1/*,
                                     -0.10316088386868619,
                                     0.13520490482776426,
                                     -0.05489235547426094,
                                     4.1434720317373253e-006,
                                     -0.00025018439997095336*/);
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
            intrinsic_params_set);

  ASSERT_EQ(0, test());
}

}
