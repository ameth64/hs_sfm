#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/synthetic/keyset_generator.hpp"

#include "hs_sfm/projective/single_camera_params_maximum_likelihood_estimator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class CorrespondencesGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
  typedef hs::sfm::synthetic::SceneGenerator<Scalar, ImageDimension>
          SceneGenerator;
  typedef typename SceneGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
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

  typedef hs::sfm::projective::SingleCameraParamsMaximumLikelihoodEstimator<
            Scalar> Estimator;
  typedef typename Estimator::Correspondence Correspondence;

  typedef EIGEN_VECTOR(Scalar, 2) Key;

public:
  typedef typename SceneGenerator::IntrinsicParams IntrinsicParams;
  typedef typename SceneGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename Estimator::CorrespondenceContainer CorrespondenceContainer;
  typedef typename Estimator::KeyCovariance KeyCovariance;

public:
  CorrespondencesGenerator(
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
      keys_generator_(image_width, image_height) {}

  Err operator() (Scalar key_stddev,
                  CorrespondenceContainer& correspondences,
                  KeyCovariance& key_covariance,
                  ExtrinsicParams& extrinsic_params_true,
                  IntrinsicParams& intrinsic_params_true) const
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

    intrinsic_params_set[0].set_skew(1e-2);
    intrinsic_params_set[0].set_principal_point_x(20);
    intrinsic_params_set[0].set_principal_point_y(30);
    intrinsic_params_set[0].set_pixel_ratio(1.0001);

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

    //生成计算P矩阵所需的对应点
    size_t number_of_tracks = tracks.size();
    key_covariance.setIdentity();
    key_covariance *= key_stddev * key_stddev;
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      if (tracks[i].size() == 1)
      {
        Key key_mean = keysets[0][tracks[i][0].second];
        Point3D point_mean = points[i];
        Correspondence correspondence(key_mean, point_mean);
        hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
          key_mean, key_covariance, correspondence.first);
        correspondences.push_back(correspondence);
      }
    }

    intrinsic_params_true = intrinsic_params_set[0];
    extrinsic_params_true = extrinsic_params_set[0];

    return 0;
  }

private:
  SceneGenerator scene_generator_;
  KeysetGenerator keys_generator_;
};

template <typename _Scalar>
class TestSingleCameraParamsMaximumLikelihoodEstimator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

private:
  typedef hs::sfm::projective::SingleCameraParamsMaximumLikelihoodEstimator<
            Scalar> Estimator;
  typedef typename Estimator::Correspondence Correspondence;

public:
  typedef typename Estimator::CorrespondenceContainer CorrespondenceContainer;
  typedef typename Estimator::ExtrinsicParams ExtrinsicParams;
  typedef typename Estimator::IntrinsicParams IntrinsicParams;
  typedef typename Estimator::KeyCovariance KeyCovariance;

public:
  Err Test(const IntrinsicParams& intrinsic_params_true,
           const ExtrinsicParams& extrinsic_params_true,
           const CorrespondenceContainer& correspondences,
           const KeyCovariance& key_covariance) const
  {
    Estimator estimator;
    Scalar focal_length_stddev = 1;
    Scalar skew_stddev = 1e-2;
    Scalar principal_point_x_stddev = 1;
    Scalar principal_point_y_stddev = 1;
    Scalar pixel_ratio_stddev = 1e-5;
    IntrinsicParams intrinsic_params_estimate;
    ExtrinsicParams extrinsic_params_estimate;
    if (estimator(correspondences, key_covariance,
                  intrinsic_params_true,
                  focal_length_stddev,
                  skew_stddev,
                  principal_point_x_stddev,
                  principal_point_y_stddev,
                  pixel_ratio_stddev,
                  intrinsic_params_estimate,
                  extrinsic_params_estimate) != 0)
    {
      return -1;
    }

    if (std::abs(intrinsic_params_true.focal_length() -
                 intrinsic_params_estimate.focal_length()) >
        focal_length_stddev * 2)
      return -1;
    if (std::abs(intrinsic_params_true.skew() -
                 intrinsic_params_estimate.skew()) >
        skew_stddev * 2)
      return -1;
    if (std::abs(intrinsic_params_true.principal_point_x() -
                 intrinsic_params_estimate.principal_point_x() >
        principal_point_x_stddev * 2))
      return -1;
    if (std::abs(intrinsic_params_true.principal_point_y() -
                 intrinsic_params_estimate.principal_point_y() >
        principal_point_y_stddev * 2))
      return -1;
    if (std::abs(intrinsic_params_true.pixel_ratio() -
                 intrinsic_params_estimate.pixel_ratio()) >
        pixel_ratio_stddev * 2)
      return -1;

    if (std::abs(extrinsic_params_true.rotation()[0] -
                 extrinsic_params_estimate.rotation()[0]) > Scalar(1e-3) ||
        std::abs(extrinsic_params_true.rotation()[1] -
                 extrinsic_params_estimate.rotation()[1]) > Scalar(1e-3) ||
        std::abs(extrinsic_params_true.rotation()[2] -
                 extrinsic_params_estimate.rotation()[2]) > Scalar(1e-3))
      return -1;
    if (!extrinsic_params_true.position().isApprox(
          extrinsic_params_estimate.position(), Scalar(0.3)))
    {
      return -1;
    }

    return 0;
  }
};

TEST(TestSingleCameraParamsMaximumLikelihoodEstimator, GeneratorTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef CorrespondencesGenerator<Scalar, ImageDimension> Generator;
  typedef TestSingleCameraParamsMaximumLikelihoodEstimator<Scalar> Tester;
  typedef Generator::CorrespondenceContainer CorrespondenceContainer;
  typedef Generator::KeyCovariance KeyCovariance;
  typedef Generator::ExtrinsicParams ExtrinsicParams;
  typedef Generator::IntrinsicParams IntrinsicParams;

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

  Generator generator(focal_length_in_metre,
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

  CorrespondenceContainer correspondences;
  KeyCovariance key_covariance;
  ExtrinsicParams extrinsic_params_true;
  IntrinsicParams intrinsic_params_true;
  Scalar key_stddev = 2;
  ASSERT_EQ(0, generator(key_stddev,
                         correspondences,
                         key_covariance,
                         extrinsic_params_true,
                         intrinsic_params_true));

  Tester tester;

  ASSERT_EQ(0, tester.Test(intrinsic_params_true, extrinsic_params_true,
                           correspondences, key_covariance));

}

//TEST(TestSingleCameraParamsMaximumLikelihoodEstimator, PriorTest)
//{
//  typedef double Scalar;
//  typedef TestSingleCameraParamsMaximumLikelihoodEstimator<Scalar> Tester;
//  typedef Tester::Correspondence Correspondence;
//  typedef Tester::CorrespondenceContainer CorrespondenceContainer;
//  typedef Tester::KeyCovariance KeyCovariance;
//  typedef Tester::ExtrinsicParams ExtrinsicParams;
//  typedef Tester::IntrinsicParams IntrinsicParams;
//
//  KeyCovariance key_covariance = KeyCovariance::Identity();
//  IntrinsicParams intrinsic_params_true(4871.7948717948721);
//  ExtrinsicParams::Rotation rotation(-0.057859973647597408,
//                                     -0.080337058970123637,
//                                     0.096320521908897089);
//  ExtrinsicParams::Position position;
//  position << -1.0056138201277076,
//              -0.38565737866881888,
//              0.0096981822538216751;
//  ExtrinsicParams extrinsic_params_true(rotation, position);
//  CorrespondenceContainer correspondences;
//  Correspondence correspondence;
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//
//  correspondence.first << -1558.8451445799994,
//                          1709.7276019096519;
//  correspondence.second << -0.59998217473680893,
//                           -0.84219614025300216,
//                           -2.4080658828714063;
//  correspondences.push_back(correspondence);
//}

}