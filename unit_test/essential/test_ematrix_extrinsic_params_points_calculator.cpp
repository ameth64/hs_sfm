#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/sfm_utility/synthetic_scene_generator.hpp"

#include "hs_sfm/essential/ematrix_5_points_calculator.hpp"
#include "hs_sfm/essential/ematrix_extrinsic_params_points_calculator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestEmatrixExtrinsicParamsPointsCalculator
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

  typedef hs::sfm::essential::EMatrix5PointsCalculator<Scalar> EMatrixCalculator;
  typedef typename EMatrixCalculator::EMatrix EMatrix;
  typedef typename EMatrixCalculator::EMatrixHypotheses EMatrixHypotheses;
  typedef typename EMatrixCalculator::HKey HKey;
  typedef typename EMatrixCalculator::HKeyPair HKeyPair;
  typedef typename EMatrixCalculator::HKeyPairContainer HKeyPairContainer;

  typedef hs::sfm::essential::EMatrixExtrinsicParamsPointsCalculator<Scalar>
          ExtrinsicParamPointsCalculator;
  typedef typename ExtrinsicParams::Rotation Rotation;
  typedef typename ExtrinsicParams::Position Position;

public:
  TestEmatrixExtrinsicParamsPointsCalculator(
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

    //生成计算E矩阵所需的匹配点
    HKeyPairContainer key_pairs;
    size_t number_of_tracks = tracks.size();
    EIGEN_MATRIX(Scalar, 2, 2) key_covariance;
    key_covariance.setIdentity();
    Scalar key_stddev = 2;
    key_covariance *= key_stddev * key_stddev;
    KMatrix K_left_inverse = intrinsic_params_set[0].GetKMatrix().inverse();
    KMatrix K_right_inverse = intrinsic_params_set[1].GetKMatrix().inverse();
    Point3DContainer points_absolute;
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      if (tracks[i].size() == 2)
      {
        size_t key_left_id = tracks[i][0].second;
        size_t key_right_id = tracks[i][1].second;
        HKeyPair key_pair;
        key_pair.first.template segment<2>(0) = keys_set[0][key_left_id];
        key_pair.second.template segment<2>(0) = keys_set[1][key_right_id];
        key_pair.first[2] = Scalar(1);
        key_pair.second[2] = Scalar(1);
        key_pair.first = K_left_inverse * key_pair.first;
        key_pair.second = K_right_inverse * key_pair.second;
        key_pairs.push_back(key_pair);
        points_absolute.push_back(points[i]);
      }
    }

    //计算由相对坐标系转换到绝对坐标系的相似变换七参数
    Rotation rotation = extrinsic_params_set[0].rotation().Inverse();
    Position translate = extrinsic_params_set[0].position();
    Scalar scale = (extrinsic_params_set[1].position() - translate).norm();

    const Rotation& rotation_absolute = extrinsic_params_set[1].rotation();
    const Position& position_absolute = extrinsic_params_set[1].position();
    Rotation rotation_relative = rotation_absolute * rotation;
    Position position_relative = (position_absolute - translate) / scale;
    position_relative = rotation.Inverse() * position_relative;
    EIGEN_MATRIX(Scalar, 3, 3) rmatrix_relative = rotation_relative;

    //计算E矩阵
    EMatrixCalculator ematrix_calculator;
    EMatrix e_matrix;
    if (ematrix_calculator(key_pairs, e_matrix) != 0)
    {
      std::cout<<"ematrix calculator failed!\n";
      return -1;
    }

    //计算外参数以及三维点
    ExtrinsicParams extrinsic_params_relative;
    Point3DContainer points_relative;
    ExtrinsicParamPointsCalculator calculator;
    if (calculator(e_matrix, key_pairs,
                   extrinsic_params_relative, points_relative) != 0)
    {
      std::cout<<"calculator failed!\n";
      return -1;
    }

    Scalar threshold = Scalar(1e-8);
    //检查相机外参数
    EIGEN_MATRIX(Scalar, 3, 3) rmatrix_relative_estimated =
      extrinsic_params_relative.rotation();
    if (!rmatrix_relative_estimated.isApprox(rmatrix_relative, threshold))
    {
      return -1;
    }
    //if (!extrinsic_params_relative.position().isApprox(position_relative,
    //                                                   threshold))
    //{
    //  return -1;
    //}

    //检查三维点
    size_t number_of_points = points_relative.size();
    for (size_t i = 0; i < number_of_points; i++)
    {
#ifdef _DEBUG
      EIGEN_MATRIX(Scalar, 3, 3) R = extrinsic_params_set[1].rotation();
      Point3D c = extrinsic_params_set[1].position();
      Point3D t = -R * c;
      EIGEN_MATRIX(Scalar, 3, 4) P;
      P.block(0, 0, 3, 3) = R;
      P.block(0, 3, 3, 1) = t;
      Point3D p = points_absolute[i];
      Scalar det = R.determinant();
      Scalar z = (P.block(2, 0, 1, 3) * p)[0] + P(2, 3);
#endif

      Point3D point_absolute = points_relative[i];
      point_absolute = rotation * point_absolute;
      point_absolute = scale * point_absolute + translate;
      if (!points_absolute[i].isApprox(point_absolute, threshold))
      {
        return -1;
      }
    }

    return 0;
  }

private:
  SceneGenerator scene_generator_;
  KeysGenerator keys_generator_;
};

TEST(TestEmatrixExtrinsicParamsPointsCalculator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestEmatrixExtrinsicParamsPointsCalculator<Scalar,
                                                     ImageDimension> Test;

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
