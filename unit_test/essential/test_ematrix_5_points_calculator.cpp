#include <iostream>

#include <gtest/gtest.h>

#include "synthetic_data_generator.hpp"

#include "hs_sfm/essential/ematrix_5_points_calculator.hpp"

namespace
{

template <typename _Scalar>
class TestEMatrix5PointsCalculator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

private:
  typedef hs::sfm::essential::EMatrix5PointsCalculator<Scalar> Calculator;
  typedef typename Calculator::EMatrix EMatrix;
  typedef typename Calculator::EMatrixHypotheses EMatrixHypotheses;
  typedef typename Calculator::HKey HKey;
  typedef typename Calculator::HKeyPair HKeyPair;
  typedef typename Calculator::HKeyPairContainer HKeyPairContainer;
public:
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef std::pair<Key, Key> KeyPair;
  typedef EIGEN_STD_VECTOR(KeyPair) KeyPairContainer;
  typedef EIGEN_MATRIX(Scalar, 3, 3) KMatrix;
private:
  typedef EIGEN_VECTOR(Scalar, 3) HLine;

public:
  TestEMatrix5PointsCalculator(){}

  Err operator() (const KeyPairContainer& key_pairs,
                  const KMatrix& K_left, const KMatrix& K_right,
                  Scalar threshold)
  {
    size_t number_of_keys = key_pairs.size();

    //生成计算E矩阵所需的匹配点
    HKeyPairContainer hkey_pairs;
    KMatrix K_left_inverse = K_left.inverse();
    KMatrix K_right_inverse = K_right.inverse();
    for (size_t i = 0; i < number_of_keys; i++)
    {
      HKeyPair hkey_pair;
      hkey_pair.first.template segment<2>(0) = key_pairs[i].first;
      hkey_pair.second.template segment<2>(0) = key_pairs[i].second;
      hkey_pair.first[2] = Scalar(1);
      hkey_pair.second[2] = Scalar(1);
      hkey_pair.first = K_left_inverse * hkey_pair.first;
      hkey_pair.second = K_right_inverse * hkey_pair.second;
      hkey_pairs.push_back(hkey_pair);
    }

    //计算E矩阵
    Calculator calculator;
    EMatrix e_matrix;
    if (calculator(hkey_pairs, e_matrix) != 0)
    {
      std::cout<<"calculator failed!\n";
      return -1;
    }

    //检查对极线误差
    Scalar mean_error = Scalar(0);
    for (size_t i = 0; i < number_of_keys; i++)
    {
      HKey hkey_left;
      hkey_left.template segment<2>(0) = key_pairs[i].first;
      hkey_left[2] = Scalar(1);
      HKey hkey_right;
      hkey_right.template segment<2>(0) = key_pairs[i].second;
      hkey_right[2] = Scalar(1);

      EIGEN_MATRIX(Scalar, 3, 3) f_matrix =
        K_right_inverse.transpose() * e_matrix * K_left_inverse;

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
    std::cout<<"mean_error:"<<mean_error<<"\n";
    if (mean_error <= threshold)
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }
};

TEST(TestEMatrix5PointsCalculator, SyntheticTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestEMatrix5PointsCalculator<Scalar> Test;

  typedef hs::sfm::essential::SyntheticDataGenerator<Scalar, ImageDimension>
          Generator;
  typedef Generator::IntrinsicParams IntrinsicParams;
  typedef Generator::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Generator::KeyPairContainer KeyPairContainer;
  typedef Generator::KMatrix KMatrix;

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

  //TODO:当生成含有镜头畸变的数据时，大量点计算的E矩阵失效。
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

  Generator generator(focal_length_in_metre,
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

  KeyPairContainer key_pairs;
  KMatrix K_left;
  KMatrix K_right;
  ASSERT_EQ(0, generator(key_pairs, K_left, K_right));

  Test test;
  ASSERT_EQ(0, test(key_pairs, K_left, K_right, Scalar(16)));
}

TEST(TestEMatrix5PointsCalculator, PriorTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestEMatrix5PointsCalculator<Scalar> Test;

  typedef Test::Key Key;
  typedef Test::KeyPair KeyPair;
  typedef Test::KeyPairContainer KeyPairContainer;
  typedef Test::KMatrix KMatrix;

  KMatrix K_left;
  K_left << 7692.3081050000001, 0, 0,
            0, 7692.3081050000001, 0,
            0, 0, 1;
  KMatrix K_right;
  K_right << 7692.3081050000001, 0, 0,
             0, 7692.3081050000001, 0,
             0, 0, 1;

  KeyPairContainer key_pairs;
  key_pairs.push_back(KeyPair(Key(3139.6198730468750, 974.15551757812500),
                              Key(3314.8464355468750, 2292.7456054687500)));
  key_pairs.push_back(KeyPair(Key(3290.5241699218750, 1903.2497558593750),
                              Key(3483.9265136718750, 3259.7797851562500)));
  key_pairs.push_back(KeyPair(Key(4760.9277343750000, 1758.7963867187500),
                              Key(4971.3750000000000, 3044.3640136718750)));
  key_pairs.push_back(KeyPair(Key(4362.1567382812500, 2254.2053222656250),
                              Key(4571.6074218750000, 3555.6372070312500)));
  key_pairs.push_back(KeyPair(Key(2854.4089355468750, 2496.6833496093750),
                              Key(3043.8037109375000, 3843.5603027343750)));
  Test test;
  ASSERT_EQ(0, test(key_pairs, K_left, K_right, Scalar(1e-8)));
}

}
