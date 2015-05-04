#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/synthetic/flight_generator.hpp"
#include "hs_sfm/synthetic/multiple_camera_keyset_generator.hpp"

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
  typedef typename KeysetGenerator::IntrinsicParams IntrinsicParams;
  typedef typename KeysetGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
private:
  typedef typename IntrinsicParams::KMatrix KMatrix;

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

    //生成计算E矩阵所需的匹配点
    HKeyPairContainer key_pairs;
    size_t number_of_tracks = tracks.size();
    EIGEN_MATRIX(Scalar, 2, 2) key_covariance;
    key_covariance.setIdentity();
    Scalar key_stddev = 2;
    key_covariance *= key_stddev * key_stddev;
    KMatrix K_left_inverse = intrinsic_params_set_[0].GetKMatrix().inverse();
    KMatrix K_right_inverse = intrinsic_params_set_[1].GetKMatrix().inverse();
    Point3DContainer points_absolute;
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      if (tracks[i].size() == 2)
      {
        size_t key_left_id = tracks[i][0].second;
        size_t key_right_id = tracks[i][1].second;
        HKeyPair key_pair;
        key_pair.first.template segment<2>(0) = keysets[0][key_left_id];
        key_pair.second.template segment<2>(0) = keysets[1][key_right_id];
        key_pair.first[2] = Scalar(1);
        key_pair.second[2] = Scalar(1);
        key_pair.first = K_left_inverse * key_pair.first;
        key_pair.second = K_right_inverse * key_pair.second;
        key_pairs.push_back(key_pair);
        points_absolute.push_back(points[i]);
      }
    }

    //计算由相对坐标系转换到绝对坐标系的相似变换七参数
    typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;
    RMatrix rotation_absolute_0 = extrinsic_params_set[0].rotation();
    RMatrix rotation = rotation_absolute_0.transpose();
    Position translate = extrinsic_params_set[0].position();
    Scalar scale = (extrinsic_params_set[1].position() - translate).norm();

    //计算真实相对外参数
    RMatrix rotation_absolute_1 = extrinsic_params_set[1].rotation();
    const Position& position_absolute = extrinsic_params_set[1].position();
    RMatrix rotation_relative = rotation_absolute_1 * rotation;
    Position position_relative = (position_absolute - translate) / scale;
    position_relative = rotation.transpose() * position_relative;

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

    ExtrinsicParams extrinsic_params_relative_sample;
    if (calculator(e_matrix, key_pairs, 20, 10,
                   extrinsic_params_relative_sample) != 0)
    {
      std::cout<<"calculator failed!\n";
      return -1;
    }

    Scalar threshold = Scalar(1e-7);
    //检查相机外参数
    EIGEN_MATRIX(Scalar, 3, 3) rmatrix_relative_estimated =
      extrinsic_params_relative.rotation();
    if (!rmatrix_relative_estimated.isApprox(rotation_relative, threshold))
    {
      std::cout<<"rotation_relative:"<<rotation_relative<<"\n";
      return -1;
    }
    if (!extrinsic_params_relative.position().isApprox(position_relative,
      threshold))
    {
      std::cout<<"position_relative:\n"<<position_relative<<"\n";
      std::cout<<"extrinsic_params_relative.position():\n"
               <<extrinsic_params_relative.position()<<"\n";
      return -1;
    }

    EIGEN_MATRIX(Scalar, 3, 3) rmatrix_relative_estimated_sample =
      extrinsic_params_relative_sample.rotation();
    if (!rmatrix_relative_estimated_sample.isApprox(
          rotation_relative, threshold))
    {
      std::cout<<"rotation_relative:"<<rotation_relative<<"\n";
      return -1;
    }
    if (!extrinsic_params_relative.position().isApprox(position_relative,
      threshold))
    {
      std::cout<<"position_relative:\n"<<position_relative<<"\n";
      std::cout<<"extrinsic_params_relative.position():\n"
               <<extrinsic_params_relative.position()<<"\n";
      return -1;
    }

    //检查三维点
    //size_t number_of_points = points_relative.size();
    //for (size_t i = 0; i < number_of_points; i++)
    //{
    //  Point3D point_absolute = points_relative[i];
    //  point_absolute = rotation * point_absolute;
    //  point_absolute = scale * point_absolute + translate;
    //  if (!points_absolute[i].isApprox(point_absolute, threshold))
    //  {
    //    std::cout<<"point_absolute:\n"<<point_absolute<<"\n";
    //    std::cout<<"points_absolute[i]:\n"<<points_absolute[i]<<"\n";
    //    return -1;
    //  }
    //}

    return 0;
  }

private:
  FlightGenerator flight_generator_;
  KeysetGenerator keys_generator_;
  IntrinsicParamsContainer intrinsic_params_set_;
};

TEST(TestEmatrixExtrinsicParamsPointsCalculator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestEmatrixExtrinsicParamsPointsCalculator<Scalar,
                                                     ImageDimension> Test;
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
