﻿#include <iostream>

#include <gtest/gtest.h>

#include "hs_test_utility/test_monte_carlo/normal_mle_simulator.hpp"

#include "hs_sfm/sfm_utility/camera_rotation_covariance.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_synthetic_data_generator.hpp"

#include "test_camera_shared_normal_mle_meta.hpp"

namespace
{

template <typename _Scalar>
class TestCameraSharedMonteCarloSimulate
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::ba::CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef
    hs::sfm::ba::CameraSharedAnalyticalJacobianMatrixCalculator<
      VectorFunction> JacobianMatrixCalculator;
  typedef typename JacobianMatrixCalculator::JacobianMatrix JacobianMatrix;
  typedef hs::test::mc::NormalMLESimulator<VectorFunction> Simulator;
  typedef typename Simulator::NoisedYVectorGenerator
                   NoisedYVectorGenerator;
  typedef typename Simulator::AnalyticXCovarianceCalculator
                   AnalyticXCovarianceCalculator;
  typedef typename Simulator::StatisticalXCovarianceCalculator
                   StatisticalXCovarianceCalculator;
  typedef typename Simulator::XVectorOptimizor
                   XVectorOptimizor;
  typedef typename Simulator::ResidualsCalculator
                   ResidualsCalculator;
  typedef typename Simulator::MahalanobisDistanceCalculator
                   MahalanobisDistanceCalculator;
  typedef typename Simulator::XVector XVector;
  typedef typename Simulator::YVector YVector;
  typedef typename Simulator::YCovarianceInverse YCovarianceInverse;
  typedef typename Simulator::XCovariance XCovariance;

  typedef EIGEN_MATRIX(Scalar, VectorFunction::params_per_key_,
                               VectorFunction::params_per_key_) KeyCovariance;

private:
  typedef hs::sfm::CameraRotaionCovarianceCalculator<Scalar>
    CameraRotationCovariance;
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic) MatrixXX;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;

public:
  TestCameraSharedMonteCarloSimulate(
    Scalar constrained_point_planar_stddev,
    Scalar constrained_point_height_stddev,
    Scalar constrained_image_rotation_stddev,
    Scalar constrained_image_position_stddev,
    Scalar constrained_camera_skew_stddev,
    Scalar constrained_camera_pixel_ratio_stddev)
    : constrained_point_planar_stddev_(constrained_point_planar_stddev),
      constrained_point_height_stddev_(constrained_point_height_stddev),
      constrained_image_rotation_stddev_(constrained_image_rotation_stddev),
      constrained_image_position_stddev_(constrained_image_position_stddev),
      constrained_camera_skew_stddev_(constrained_camera_skew_stddev),
      constrained_camera_pixel_ratio_stddev_(
        constrained_camera_pixel_ratio_stddev) {}

  Err Test(const VectorFunction& vector_function,
           const XVector& true_x,
           const KeyCovariance& key_covariance) const
  {
    YVector true_y;
    if (vector_function(true_x, true_y) != 0)
    {
      return -1;
    }

    XVector near_x;
    if (GetNearX(vector_function, true_x, near_x) != 0)
    {
      return -1;
    }

    XVectorOptimizor x_optimizor(near_x, 100,
                                 Scalar(1e-6),
                                 Scalar(1e-10),
                                 Scalar(1e-10));

    NoisedYVectorGenerator noised_y_generator;

    AnalyticXCovarianceCalculator analytic_x_covariance_calculator;

    ResidualsCalculator residuals_calculator;

    MahalanobisDistanceCalculator mahalanobis_distance_calculator;

    YCovarianceInverse y_covariance_inverse;
    if (GenerateYCovarianceInverse(vector_function, key_covariance,
                                   y_covariance_inverse) != 0)
    {
      return -1;
    }

    StatisticalXCovarianceCalculator statistical_x_covariance_calculator;
    XCovariance analytic_x_covariance;
    XCovariance statistical_x_covariance;
    Scalar residual_mean;

    Simulator simulator;
    if (simulator(vector_function,
                  x_optimizor,
                  residuals_calculator,
                  mahalanobis_distance_calculator,
                  noised_y_generator,
                  analytic_x_covariance_calculator,
                  true_y,
                  y_covariance_inverse,
                  true_x,
                  5,
                  statistical_x_covariance_calculator,
                  residual_mean,
                  analytic_x_covariance,
                  statistical_x_covariance) != 0)
    {
      return -1;
    }
    residual_mean = std::sqrt(residual_mean / Scalar(true_y.rows()));

    //计算本质参数数
    JacobianMatrixCalculator jacobian;
    JacobianMatrix jacobian_matrix;
    if (jacobian(vector_function, true_x, jacobian_matrix) != 0)
    {
      std::cout<<"jacobian failed!\n";
      return -1;
    }
    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    MatrixXX dense_jacobian_matrix(y_size, x_size);
    for (Index i = 0; i < y_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        dense_jacobian_matrix(i, j) = jacobian_matrix.coeff(i, j);
      }
    }
    Eigen::FullPivLU<MatrixXX> lu(dense_jacobian_matrix);
    Index rank = lu.rank();

    Scalar expected_residual_mean =
      std::sqrt(Scalar(1) - Scalar(rank) / Scalar(y_size));

    Err result = 0;
    if (std::abs(residual_mean - expected_residual_mean) > 5e-2)
    {
      std::cout<<"residual mean too large\n";
      std::cout<<"expected residual mean:"<<expected_residual_mean<<"\n";
      std::cout<<"residual mean:"<<residual_mean<<"\n";
      result = -1;
    }

    return result;
  }

private:
  Err GenerateYCovarianceInverse(
    const VectorFunction& vector_function,
    const KeyCovariance& key_covariance,
    YCovarianceInverse& y_covariance_inverse) const
  {
    size_t number_of_keys = size_t(vector_function.number_of_keys());
    y_covariance_inverse.Clear();
    y_covariance_inverse.SetKeysUniformCovariance(key_covariance,
                                                  number_of_keys);
    auto itr_point_constriant = vector_function.point_constraints().begin();
    auto itr_point_constriant_end = vector_function.point_constraints().end();
    for (; itr_point_constriant != itr_point_constriant_end;
         ++itr_point_constriant)
    {
      MatrixXX point_block;
      Scalar planar_constraint =
        Scalar(1) / (constrained_point_planar_stddev_ *
                     constrained_point_planar_stddev_);
      Scalar height_constraint =
        Scalar(1) / (constrained_point_height_stddev_ *
                     constrained_point_height_stddev_);
      if (itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_X] &&
          itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Y] &&
          !itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Z])
      {
        y_covariance_inverse.AddConstraint(planar_constraint);
        y_covariance_inverse.AddConstraint(planar_constraint);
      }
      else if (itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_X] &&
               itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Y] &&
               itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Z])
      {
        y_covariance_inverse.AddConstraint(planar_constraint);
        y_covariance_inverse.AddConstraint(planar_constraint);
        y_covariance_inverse.AddConstraint(height_constraint);
      }
    }

    auto itr_image_constraint = vector_function.image_constraints().begin();
    auto itr_image_constraint_end = vector_function.image_constraints().end();
    const Scalar pi = Scalar(3.141592653);
    for (; itr_image_constraint != itr_image_constraint_end;
         ++itr_image_constraint)
    {
      MatrixXX image_block;
      if (itr_image_constraint->mask[hs::sfm::ba::IMAGE_CONSTRAIN_ROTATION] &&
          itr_image_constraint->mask[hs::sfm::ba::IMAGE_CONSTRAIN_POSITION_X] &&
          itr_image_constraint->mask[hs::sfm::ba::IMAGE_CONSTRAIN_POSITION_Y] &&
          itr_image_constraint->mask[hs::sfm::ba::IMAGE_CONSTRAIN_POSITION_Z])
      {
        Scalar rotation_stddev_radian =
          constrained_image_rotation_stddev_ / Scalar(180) * pi;
        Scalar rotation_constraint =
          Scalar(1) / (rotation_stddev_radian * rotation_stddev_radian);
        Scalar position_constraint =
          Scalar(1) / (constrained_image_position_stddev_ *
                       constrained_image_position_stddev_);
        y_covariance_inverse.AddConstraint(rotation_constraint);
        y_covariance_inverse.AddConstraint(rotation_constraint);
        y_covariance_inverse.AddConstraint(rotation_constraint);
        y_covariance_inverse.AddConstraint(position_constraint);
        y_covariance_inverse.AddConstraint(position_constraint);
        y_covariance_inverse.AddConstraint(position_constraint);
      }
    }

    auto itr_camera_constraint = vector_function.camera_constraints().begin();
    auto itr_camera_constraint_end = vector_function.camera_constraints().end();
    for (; itr_camera_constraint != itr_camera_constraint_end;
         ++itr_camera_constraint)
    {
      if (itr_camera_constraint->radial_mask.none() &&
          itr_camera_constraint->decentering_mask.none() &&
          !itr_camera_constraint->intrinsic_mask[
            hs::sfm::ba::INTRINSIC_CONSTRAIN_FOCAL_LENGTH] &&
          itr_camera_constraint->intrinsic_mask[
            hs::sfm::ba::INTRINSIC_CONSTRAIN_SKEW] &&
          !itr_camera_constraint->intrinsic_mask[
            hs::sfm::ba::INTRINSIC_CONSTRAIN_PRINCIPAL_X] &&
          !itr_camera_constraint->intrinsic_mask[
            hs::sfm::ba::INTRINSIC_CONSTRAIN_PRINCIPAL_Y] &&
          itr_camera_constraint->intrinsic_mask[
            hs::sfm::ba::INTRINSIC_CONSTRAIN_PIXEL_RATIO])
      {
        Scalar skew_constraint =
          Scalar(1) / (constrained_camera_skew_stddev_ *
                       constrained_camera_skew_stddev_);
        Scalar pixel_ratio_constaint =
          Scalar(1) / (constrained_camera_pixel_ratio_stddev_ *
                       constrained_camera_pixel_ratio_stddev_);
        y_covariance_inverse.AddConstraint(skew_constraint);
        y_covariance_inverse.AddConstraint(pixel_ratio_constaint);
      }
    }
    return 0;
  }

  Err GetNearX(const VectorFunction& vector_function,
               const XVector& true_x,
               XVector& near_x) const
  {
    Index x_size = true_x.rows();
    near_x.resize(x_size);

    Index number_of_points = vector_function.number_of_points();
    Scalar point_planar_stddev = Scalar(1);
    Scalar point_height_stddev = Scalar(3);
    for (Index i = 0; i < number_of_points; i++)
    {
      Vector3 point = true_x.segment(i * VectorFunction::params_per_point_,
                                     VectorFunction::params_per_point_);
      Matrix33 point_covariance = Matrix33::Zero();
      point_covariance(0, 0) = point_planar_stddev * point_planar_stddev;
      point_covariance(1, 1) = point_planar_stddev * point_planar_stddev;
      point_covariance(2, 2) = point_height_stddev * point_height_stddev;
      Vector3 noised_point;
      hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
        point, point_covariance, noised_point);
      near_x.segment(i * VectorFunction::params_per_point_,
                     VectorFunction::params_per_point_) = noised_point;
    }

    Index number_of_images = vector_function.number_of_images();
    Scalar image_rotation_stddev = Scalar(1);
    Scalar image_position_stddev = Scalar(5);
    Index point_params_size = vector_function.GetPointParamsSize();
    CameraRotationCovariance rotation_covariance_calculator;
    for (Index i = 0; i < number_of_images; i++)
    {
      Vector3 camera_rotation =
        true_x.segment(point_params_size +
                       i * VectorFunction::extrinsic_params_per_image_, 3);
      Matrix33 camera_rotation_covariance;
      if (rotation_covariance_calculator(camera_rotation,
                                         image_rotation_stddev,
                                         image_rotation_stddev,
                                         image_rotation_stddev,
                                         camera_rotation_covariance) != 0)
      {
        return -1;
      }
      Vector3 noised_camera_rotation;
      hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
        camera_rotation, camera_rotation_covariance, noised_camera_rotation);
      near_x.segment(point_params_size +
                     i * VectorFunction::extrinsic_params_per_image_, 3) =
        noised_camera_rotation;

      Vector3 camera_pos =
        true_x.segment(point_params_size +
                       i * VectorFunction::extrinsic_params_per_image_ + 3, 3);
      Matrix33 camera_pos_covariance = Matrix33::Zero();
      camera_pos_covariance(0, 0) =
        image_position_stddev * image_position_stddev;
      camera_pos_covariance(1, 1) =
        image_position_stddev * image_position_stddev;
      camera_pos_covariance(2, 2) =
        image_position_stddev * image_position_stddev;
      Vector3 noised_camera_pos;
      hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
        camera_pos, camera_pos_covariance, noised_camera_pos);
      near_x.segment(point_params_size +
                     i * VectorFunction::extrinsic_params_per_image_ + 3, 3) =
        noised_camera_pos;
    }

    Index number_of_cameras = vector_function.number_of_cameras();
    Index extrinsic_params_size = vector_function.GetExtrinsicParamsSize();
    Index intrinsic_params_size_per_camera =
      vector_function.GetIntrinsicParamsSizePerCamera();
    const hs::sfm::ba::IntrinsicComputationsMask& intrinsic_computations_mask =
      vector_function.intrinsic_computations_mask();
    Scalar focal_length_stddev = 10;
    for (Index i = 0; i < number_of_cameras; i++)
    {
      Index offset = point_params_size + extrinsic_params_size +
                     i * intrinsic_params_size_per_camera;
      if (intrinsic_computations_mask[
            hs::sfm::ba::COMPUTE_RADIAL_DISTORTION])
      {
        near_x[offset + 0] = Scalar(0);
        near_x[offset + 1] = Scalar(0);
        near_x[offset + 2] = Scalar(0);
        offset += 3;
      }
      if (intrinsic_computations_mask[
            hs::sfm::ba::COMPUTE_DECENTERING_DISTORTION])
      {
        near_x[offset + 0] = Scalar(0);
        near_x[offset + 1] = Scalar(0);
        offset += 2;
      }
      if (intrinsic_computations_mask[
            hs::sfm::ba::COMPUTE_INTRINSIC_PARAMS])
      {
        hs::math::random::NormalRandomVar<Scalar, 1>::Generate(
          true_x[offset], focal_length_stddev, near_x[offset]);
        near_x[offset + 1] = Scalar(0);
        near_x[offset + 2] = Scalar(0);
        near_x[offset + 3] = Scalar(0);
        near_x[offset + 4] = Scalar(1);
      }
    }

    return 0;
  }

private:
  Scalar constrained_point_planar_stddev_;
  Scalar constrained_point_height_stddev_;
  Scalar constrained_image_rotation_stddev_;
  Scalar constrained_image_position_stddev_;
  Scalar constrained_camera_skew_stddev_;
  Scalar constrained_camera_pixel_ratio_stddev_;
};

TEST(TestCameraSharedMonteCarloSimulate, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestCameraSharedMonteCarloSimulate<Scalar>
          Tester;
  typedef Tester::KeyCovariance KeyCovariance;

  typedef hs::sfm::ba::CameraSharedSyntheticDataGenerator<Scalar,
                                                          ImageDimension>
          DataGenerator;
  typedef DataGenerator::FlightGenerator FlightGenerator;
  typedef DataGenerator::FlightGeneratorContainer FlightGeneratorContainer;
  typedef DataGenerator::IntrinsicParams IntrinsicParams;
  typedef DataGenerator::IntrinsicParamsContainer IntrinsicParamsContainer;

  typedef hs::sfm::ba::CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef VectorFunction::Index Index;
  typedef VectorFunction::XVector XVector;
  typedef VectorFunction::YVector YVector;

  FlightGeneratorContainer flight_generators;
  IntrinsicParamsContainer intrinsic_params_set;

  Scalar focal_length_in_metre_0 = 0.018858358970276164;
  size_t number_of_strips_0 = 3;
  size_t number_of_cameras_in_strip_0 = 3;
  Scalar ground_resolution_0 = 0.1;
  ImageDimension image_width_0 = 6000;
  ImageDimension image_height_0 = 4000;
  Scalar pixel_size_0 = 0.0000039;
  size_t number_of_points_0 = 2000;
  Scalar lateral_overlap_ratio_0 = 0.6;
  Scalar longitudinal_overlap_ratio_0 = 0.8;
  Scalar scene_max_height_0 = 50;
  Scalar camera_height_stddev_0 = 2;
  Scalar camera_planar_stddev_0 = 2;
  Scalar camera_rotation_stddev_0 = 10;
  FlightGenerator flight_generator_0(
    focal_length_in_metre_0,
    number_of_strips_0,
    number_of_cameras_in_strip_0,
    ground_resolution_0,
    image_width_0,
    image_height_0,
    pixel_size_0,
    number_of_points_0,
    lateral_overlap_ratio_0,
    longitudinal_overlap_ratio_0,
    scene_max_height_0,
    camera_height_stddev_0,
    camera_planar_stddev_0,
    camera_rotation_stddev_0);
  flight_generators.push_back(flight_generator_0);
  IntrinsicParams intrinsic_params_0(focal_length_in_metre_0 / pixel_size_0,
                                     0,
                                     -42.4095312016,
                                     -31.699212823,
                                     1,
                                     -0.0050490462006048831,
                                     0.031293804298609881,
                                     -0.030794820960442223,
                                     -0.00055376548320189127,
                                     -0.00049877717768381476);
  intrinsic_params_set.push_back(intrinsic_params_0);

  Scalar focal_length_in_metre_1 = 0.02995452167701055;
  size_t number_of_strips_1 = 4;
  size_t number_of_cameras_in_strip_1 = 3;
  Scalar ground_resolution_1 = 0.1;
  ImageDimension image_width_1 = 6000;
  ImageDimension image_height_1 = 4000;
  Scalar pixel_size_1 = 0.0000039;
  size_t number_of_points_1 = 2000;
  Scalar lateral_overlap_ratio_1 = 0.6;
  Scalar longitudinal_overlap_ratio_1 = 0.8;
  Scalar scene_max_height_1 = 50;
  Scalar camera_height_stddev_1 = 2;
  Scalar camera_planar_stddev_1 = 2;
  Scalar camera_rotation_stddev_1 = 10;
  FlightGenerator flight_generator_1(
    focal_length_in_metre_1,
    number_of_strips_1,
    number_of_cameras_in_strip_1,
    ground_resolution_1,
    image_width_1,
    image_height_1,
    pixel_size_1,
    number_of_points_1,
    lateral_overlap_ratio_1,
    longitudinal_overlap_ratio_1,
    scene_max_height_1,
    camera_height_stddev_1,
    camera_planar_stddev_1,
    camera_rotation_stddev_1);
  flight_generators.push_back(flight_generator_1);
  IntrinsicParams intrinsic_params_1(focal_length_in_metre_1 / pixel_size_1,
                                     0,
                                     -21.669436058,
                                     -44.8644764322,
                                     1,
                                     -0.02529179096221609,
                                     0.23762413973445157,
                                     -0.64208397668697237,
                                     -0.0020605099808780948,
                                     -0.00028706423764766859);
  intrinsic_params_set.push_back(intrinsic_params_1);

  Scalar focal_length_in_metre_2 = 0.019056097774998712;
  size_t number_of_strips_2 = 2;
  size_t number_of_cameras_in_strip_2 = 5;
  Scalar ground_resolution_2 = 0.1;
  ImageDimension image_width_2 = 6000;
  ImageDimension image_height_2 = 4000;
  Scalar pixel_size_2 = 0.0000039;
  size_t number_of_points_2 = 2000;
  Scalar lateral_overlap_ratio_2 = 0.6;
  Scalar longitudinal_overlap_ratio_2 = 0.8;
  Scalar scene_max_height_2 = 50;
  Scalar camera_height_stddev_2 = 2;
  Scalar camera_planar_stddev_2 = 2;
  Scalar camera_rotation_stddev_2 = 10;
  FlightGenerator flight_generator_2(
    focal_length_in_metre_2,
    number_of_strips_2,
    number_of_cameras_in_strip_2,
    ground_resolution_2,
    image_width_2,
    image_height_2,
    pixel_size_2,
    number_of_points_2,
    lateral_overlap_ratio_2,
    longitudinal_overlap_ratio_2,
    scene_max_height_2,
    camera_height_stddev_2,
    camera_planar_stddev_2,
    camera_rotation_stddev_2);
  flight_generators.push_back(flight_generator_2);
  IntrinsicParams intrinsic_params_2(focal_length_in_metre_2 / pixel_size_2,
                                     0,
                                     -35.2052431556,
                                     -16.4262220759,
                                     1,
                                     -0.10316088386868619,
                                     0.13520490482776426,
                                     -0.05489235547426094,
                                     4.1434720317373253e-006,
                                     -0.00025018439997095336);
  intrinsic_params_set.push_back(intrinsic_params_2);

  Scalar flight_longitudinal_overlap_ratio = 0.8;
  Scalar flight_lateral_overlap_ratio = 0.2;
  Scalar north_west_angle = 60.0;
  Scalar north_west_angle_stddev = 10.0;
  Scalar offset_stddev = 15.0;
  size_t number_of_points = 100;
  size_t number_of_planar_constrained_points = 4;
  size_t number_of_full_constrained_points = 6;
  size_t number_of_constrained_images = 10;
  size_t number_of_constrained_cameras = 2;

  Scalar constrianed_point_planar_stddev = 0.05;
  Scalar constrianed_point_height_stddev = 0.1;
  Scalar constrianed_image_rotation_stddev = 0.5;
  Scalar constrianed_image_position_stddev = 1;
  Scalar constrianed_camera_skew_stddev = 0.001;
  Scalar constrianed_camera_pixel_ratio_stddev = 0.001;
  Tester tester(constrianed_point_planar_stddev,
                constrianed_point_height_stddev,
                constrianed_image_rotation_stddev,
                constrianed_image_position_stddev,
                constrianed_camera_skew_stddev,
                constrianed_camera_pixel_ratio_stddev);
  KeyCovariance feature_covariance = KeyCovariance::Identity();
  //feature_covariance *= Scalar(4);

  int number_of_fix_masks = (1 << hs::sfm::ba::NUMBER_OF_FIX_PARAMS) - 1;
  for (int i = 0; i < number_of_fix_masks; i++)
  {
    hs::sfm::ba::FixMask fix_mask(i);
    DataGenerator data_generator(flight_longitudinal_overlap_ratio,
                                 flight_lateral_overlap_ratio,
                                 north_west_angle,
                                 north_west_angle_stddev,
                                 offset_stddev,
                                 flight_generators,
                                 number_of_points,
                                 number_of_planar_constrained_points,
                                 number_of_full_constrained_points,
                                 number_of_constrained_images,
                                 number_of_constrained_cameras,
                                 intrinsic_params_set,
                                 fix_mask);

    VectorFunction vector_function;
    XVector x;
    YVector y_synthetic, y;

    ASSERT_EQ(0, data_generator(vector_function, x, y_synthetic));
    ASSERT_EQ(0, vector_function(x, y));

    ASSERT_EQ(0, tester.Test(vector_function, x, feature_covariance));
  }
}

}
