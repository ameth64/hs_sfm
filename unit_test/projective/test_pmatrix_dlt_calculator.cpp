#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/synthetic/flight_generator.hpp"
#include "hs_sfm/synthetic/multiple_camera_keyset_generator.hpp"

#include "hs_sfm/projective/pmatrix_dlt_calculator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestPMatrixDLTCalculator
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
  typedef typename KeysetGenerator::IntrinsicParams IntrinsicParams;
  typedef typename KeysetGenerator::IntrinsicParamsContainer
          IntrinsicParamsContainer;

  typedef hs::sfm::projective::PMatrixDLTCalculator<Scalar> Calculator;
  typedef typename Calculator::Key Key;
  typedef typename Calculator::Point Point;
  typedef typename Calculator::Correspondence Correspondence;
  typedef typename Calculator::CorrespondenceContainer CorrespondenceContainer;
  typedef typename Calculator::PMatrix PMatrix;
  typedef EIGEN_VECTOR(Scalar, 3) HKey;
  typedef EIGEN_VECTOR(Scalar, 4) HPoint;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;

public:
  TestPMatrixDLTCalculator(
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
    Scalar north_west_angle)
    : flight_generator_(focal_length_in_metre,
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
                       camera_rot_stdddev) {}

  Err operator() ()
  {
    //生成相机参数和三维点
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    Point3DContainer points;
    if (flight_generator_(extrinsic_params_set,
                          images, points) != 0)
    {
      std::cout<<"scene generator failed!\n";
      return -1;
    }

    IntrinsicParamsContainer intrinsic_params_set(
      extrinsic_params_set.size(), flight_generator_.GetFocalLengthInPixel());
    std::vector<size_t> image_intrinsic_map;
    image_intrinsic_map.push_back(0);

    //生成特征点
    KeysetContainer keysets;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;
    if (keys_generator_(intrinsic_params_set,
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

    //生成计算P矩阵所需的对应点
    CorrespondenceContainer correspondences;
    size_t number_of_tracks = tracks.size();
    Scalar key_stddev = 2;
    EIGEN_MATRIX(Scalar, 2, 2) key_covariance;
    key_covariance.setIdentity();
    key_covariance *= key_stddev * key_stddev;
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      if (tracks[i].size() == 1)
      {
        Key key_mean = keysets[0][tracks[i][0].second];
        Point point_mean = points[i];
        Correspondence correspondence(key_mean, point_mean);
        hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
          key_mean, key_covariance, correspondence.first);
        correspondences.push_back(correspondence);
      }
    }

    //计算P矩阵
    Calculator calculator;
    PMatrix p_matrix;
    if (calculator(correspondences, p_matrix) != 0)
    {
      std::cout<<"calculator failed!\n";
      return -1;
    }

    //检查重投影误差
    Scalar mean_error = Scalar(0);
    size_t number_of_correspondences = correspondences.size();
    for (size_t i = 0; i < number_of_correspondences; i++)
    {
      HPoint hpoint;
      hpoint.template segment<3>(0) = correspondences[i].second;
      hpoint[3] = Scalar(1);
      HKey hkey = p_matrix * hpoint;
      hkey /= hkey(2);
      Scalar error = (hkey.template segment<2>(0) -
                      correspondences[i].first).norm();
      mean_error += error;
    }
    mean_error /= Scalar(number_of_correspondences);

    if (mean_error <= key_stddev + 1)
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }

private:
  FlightGenerator flight_generator_;
  KeysetGenerator keys_generator_;
};

TEST(TestPMatrixDLTCalculator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestPMatrixDLTCalculator<Scalar, ImageDimension> Test;

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
            north_west_angle);

  ASSERT_EQ(0, test());

}

}
