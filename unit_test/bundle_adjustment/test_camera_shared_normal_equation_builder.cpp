﻿#include <iostream>

#include <gtest/gtest.h>

#define IMAGE_TEST
#ifdef IMAGE_TEST
#include "hs_image_io/whole_io/image_data.hpp"
#include "hs_image_io/whole_io/image_io.hpp"
#include <fstream>
#endif

#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_synthetic_data_generator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_noised_y_generator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_normal_equation_builder.hpp"

namespace
{

template <typename _Scalar>
class TestCameraSharedNormalEquationBuilder
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef hs::sfm::ba::CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;

  typedef hs::sfm::ba::CameraSharedNormalEquationBuilder<Scalar>
          NormalEquationBuilder;
  typedef typename NormalEquationBuilder::YCovarianceInverse
                   YCovarianceInverse;
  typedef typename NormalEquationBuilder::NormalMatrix NormalMatrix;
  typedef typename NormalEquationBuilder::Gradient Gradient;
  typedef typename NormalEquationBuilder::Residuals Residuals;

  typedef hs::sfm::ba::CameraSharedAnalyticalJacobianMatrixCalculator<
            VectorFunction> JacobianMatrixCalculator;
  typedef typename JacobianMatrixCalculator::JacobianMatrix JacobianMatrix;

  typedef hs::sfm::ba::CameraSharedNoisedYGenerator<Scalar>
          NoisedYGenerator;

  typedef EIGEN_MATRIX(Scalar, VectorFunction::params_per_key_,
                               VectorFunction::params_per_key_)
          KeyCovariance;

private:
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic) MatrixXX;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) VectorX;

public:
  TestCameraSharedNormalEquationBuilder(
    Scalar point_planar_stddev,
    Scalar point_height_stddev,
    Scalar image_rotation_stddev,
    Scalar image_position_stddev,
    Scalar camera_skew_stddev,
    Scalar camera_pixel_ratio_stddev)
    : point_planar_stddev_(point_planar_stddev),
      point_height_stddev_(point_height_stddev),
      image_rotation_stddev_(image_rotation_stddev),
      image_position_stddev_(image_position_stddev),
      camera_skew_stddev_(camera_skew_stddev),
      camera_pixel_ratio_stddev_(camera_pixel_ratio_stddev) {}

  Err Test(const VectorFunction& vector_function,
           const XVector& x,
           const KeyCovariance& key_covariance)
  {
    YCovarianceInverse y_covariance_inverse;
    if (GenerateYCovarianceInverse(vector_function,
                                   key_covariance,
                                   y_covariance_inverse) != 0)
    {
      return -1;
    }

    MatrixXX y_covariance_inverse_dense;
    if (GenerateDenseYCovarianceInverse(vector_function,
                                        key_covariance,
                                        y_covariance_inverse_dense) != 0)
    {
      return -1;
    }

    if (CompareYCovarianceInverse(y_covariance_inverse,
                                  y_covariance_inverse_dense) != 0)
    {
      return -1;
    }

    Residuals residuals;
    if (GenerateResiduals(vector_function,
                          x,
                          y_covariance_inverse,
                          residuals) != 0)
    {
      return -1;
    }

    JacobianMatrix jacobian_matrix;
    if (CalculateJacobianMatrix(vector_function, x, jacobian_matrix) != 0)
    {
      return -1;
    }

    MatrixXX jacobian_matrix_dense;
    if (GenerateDenseJacobianMatrix(jacobian_matrix,
                                    jacobian_matrix_dense) != 0)
    {
      return -1;
    }

    MatrixXX normal_matrix_dense;
    VectorX gradient_dense;
    if (GenerateDenseNormalEquation(jacobian_matrix_dense,
                                    residuals,
                                    y_covariance_inverse_dense,
                                    normal_matrix_dense,
                                    gradient_dense) != 0)
    {
      return -1;
    }

    NormalMatrix normal_matrix;
    Gradient gradient;
    if (CalculateNormalEquation(jacobian_matrix,
                                residuals,
                                y_covariance_inverse,
                                normal_matrix,
                                gradient) != 0)
    {
      return -1;
    }

    if (CompareNormalEquation(normal_matrix,
                              gradient,
                              normal_matrix_dense,
                              gradient_dense) != 0)
    {
      return -1;
    }

    return 0;
  }

private:
  Err GenerateYCovarianceInverse(const VectorFunction& vector_function,
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
        Scalar(1) / (point_planar_stddev_ * point_planar_stddev_);
      Scalar height_constraint =
        Scalar(1) / (point_height_stddev_ * point_height_stddev_);
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
          image_rotation_stddev_ / Scalar(180) * pi;
        Scalar rotation_constraint =
          Scalar(1) / (rotation_stddev_radian * rotation_stddev_radian);
        Scalar position_constraint =
          Scalar(1) / (image_position_stddev_ * image_position_stddev_);
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
          Scalar(1) / (camera_skew_stddev_ * camera_skew_stddev_);
        Scalar pixel_ratio_constaint =
          Scalar(1) / (camera_pixel_ratio_stddev_ * camera_pixel_ratio_stddev_);
        y_covariance_inverse.AddConstraint(skew_constraint);
        y_covariance_inverse.AddConstraint(pixel_ratio_constaint);
      }
    }

    return 0;
  }

  Err GenerateDenseYCovarianceInverse(
    const VectorFunction& vector_function,
    const KeyCovariance& key_covariance,
    MatrixXX& y_covariance_inverse_dense) const
  {
    Index y_size = vector_function.GetYSize();
    y_covariance_inverse_dense.resize(y_size, y_size);
    y_covariance_inverse_dense.setIdentity();
    Index number_of_keys = vector_function.number_of_keys();
    MatrixXX key_covariance_inverse = key_covariance.inverse();
    for (Index i = 0; i < number_of_keys; i++)
    {
      y_covariance_inverse_dense.block(i * VectorFunction::params_per_key_,
                                       i * VectorFunction::params_per_key_,
                                       VectorFunction::params_per_key_,
                                       VectorFunction::params_per_key_) =
        key_covariance_inverse;
    }

    Index offset = vector_function.GetYKeysSize();
    auto itr_point_constriant = vector_function.point_constraints().begin();
    auto itr_point_constriant_end = vector_function.point_constraints().end();
    for (; itr_point_constriant != itr_point_constriant_end;
         ++itr_point_constriant)
    {
      MatrixXX point_block;
      if (itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_X] &&
          itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Y] &&
          !itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Z])
      {
        point_block.resize(2, 2);
        point_block.setIdentity();
        point_block /= point_planar_stddev_ * point_planar_stddev_;
        y_covariance_inverse_dense.block(offset, offset, 2, 2) = point_block;
        offset += 2;
      }
      else if (itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_X] &&
               itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Y] &&
               itr_point_constriant->mask[hs::sfm::ba::POINT_CONSTRAIN_Z])
      {
        point_block.resize(3, 3);
        point_block.setIdentity();
        point_block(0, 0) /= point_planar_stddev_ * point_planar_stddev_;
        point_block(1, 1) /= point_planar_stddev_ * point_planar_stddev_;
        point_block(2, 2) /= point_height_stddev_ * point_height_stddev_;
        y_covariance_inverse_dense.block(offset, offset, 3, 3) = point_block;
        offset += 3;
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
        image_block.resize(6, 6);
        image_block.setIdentity();
        Scalar rotation_stddev_radian =
          image_rotation_stddev_ / Scalar(180) * pi;
        image_block.block(0, 0, 3, 3) /=
          rotation_stddev_radian * rotation_stddev_radian;
        image_block.block(3, 3, 3, 3) /=
          image_position_stddev_ * image_position_stddev_;
        y_covariance_inverse_dense.block(offset, offset, 6, 6) = image_block;
        offset += 6;
      }
    }

    auto itr_camera_constraint = vector_function.camera_constraints().begin();
    auto itr_camera_constraint_end = vector_function.camera_constraints().end();
    for (; itr_camera_constraint != itr_camera_constraint_end;
         ++itr_camera_constraint)
    {
      MatrixXX camera_block;
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
        camera_block.resize(2, 2);
        camera_block.setIdentity();
        camera_block(0, 0) /=
          camera_skew_stddev_ * camera_skew_stddev_;
        camera_block(1, 1) /=
          camera_pixel_ratio_stddev_ * camera_pixel_ratio_stddev_;
        y_covariance_inverse_dense.block(offset, offset, 2, 2) = camera_block;
        offset += 2;
      }
    }

    return 0;
  }

  Err CompareYCovarianceInverse(
    const YCovarianceInverse& y_covariance_inverse,
    const MatrixXX& y_covariance_inverse_dense) const
  {
    Index x_size = y_covariance_inverse_dense.cols();
    Index y_size = y_covariance_inverse_dense.rows();
    const Scalar threshold = Scalar(1e-8);

    Err result = 0;
    for (Index i = 0; i < y_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        Scalar value = y_covariance_inverse.coeff(i, j);
        Scalar value_dense = y_covariance_inverse_dense(i, j);
        Scalar error = std::abs(value - value_dense);
        if (error > threshold)
        {
          std::cout<<"error:"<<error<<"\n";
          std::cout<<"      y_covariance_inverse["<<i<<","<<j<<"]="
                   <<value<<"\n";
          std::cout<<"y_covariance_inverse_dense["<<i<<","<<j<<"]="
                   <<value_dense<<"\n";
          result = -1;
        }
      }
    }
    return result;
  }

  Err GenerateResiduals(const VectorFunction& vector_function,
                        const XVector& x,
                        const YCovarianceInverse& y_covariance_inverse,
                        Residuals& residuals) const
  {
    YVector y;
    if (vector_function(x, y) != 0) return -1;
    NoisedYGenerator generator;
    YVector noised_y;
    if (generator(y, y_covariance_inverse, noised_y) != 0) return -1;
    residuals = y - noised_y;
    return 0;
  }

  Err CalculateJacobianMatrix(const VectorFunction& vector_function,
                              const XVector& x,
                              JacobianMatrix& jacobian_matrix) const
  {
    JacobianMatrixCalculator calculator;
    return calculator(vector_function, x, jacobian_matrix);
  }

  Err GenerateDenseJacobianMatrix(const JacobianMatrix& jacobian_matrix,
                                  MatrixXX& jacobian_matrix_dense) const
  {
    Index x_size = jacobian_matrix.GetXSize();
    Index y_size = jacobian_matrix.GetYSize();
    jacobian_matrix_dense.resize(y_size, x_size);
    for (Index i = 0; i < y_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        jacobian_matrix_dense(i, j) = jacobian_matrix.coeff(i, j);
      }
    }
    return 0;
  }

  Err CalculateNormalEquation(const JacobianMatrix& jacobian_matrix,
                              const Residuals& residuals,
                              const YCovarianceInverse& y_covariance_inverse,
                              NormalMatrix& normal_matrix,
                              Gradient& gradient) const
  {
    NormalEquationBuilder builder;
    return builder(jacobian_matrix,
                   residuals,
                   y_covariance_inverse,
                   normal_matrix,
                   gradient);
  }

  Err GenerateDenseNormalEquation(const MatrixXX& jacobian_matrix_dense,
                                  const VectorX& residuals,
                                  const MatrixXX& y_covariance_inverse_dense,
                                  MatrixXX& normal_matrix_dense,
                                  VectorX& gradient_dense) const
  {
    Index x_size = jacobian_matrix_dense.cols();
    Index y_size = jacobian_matrix_dense.rows();
    if (y_size != residuals.rows()) return -1;
    if (y_size != y_covariance_inverse_dense.rows() ||
        y_size != y_covariance_inverse_dense.cols()) return -1;
    normal_matrix_dense.resize(x_size, x_size);
    gradient_dense.resize(x_size);
    normal_matrix_dense = jacobian_matrix_dense.transpose() *
                          y_covariance_inverse_dense *
                          jacobian_matrix_dense;
    gradient_dense = jacobian_matrix_dense.transpose() *
                     y_covariance_inverse_dense *
                     residuals;
#ifdef IMAGE_TEST
    SaveMatrixImage(normal_matrix_dense, "normal_matrix_dense.png");
    SaveMatrixImage(jacobian_matrix_dense, "jacobian_matrix_dense.png");
#endif
    return 0;
  }

#ifdef IMAGE_TEST
  Err SaveMatrixImage(const MatrixXX& matrix,
                      const std::string& image_path) const
  {
    typedef hs::imgio::whole::ImageData::Byte Byte;
    hs::imgio::whole::ImageData image_data;
    Index x_size = matrix.cols();
    Index y_size = matrix.rows();
    image_data.CreateImage(int(x_size), int(y_size), 1);
    image_data.SetColorType(IMAGE_GRAYSCALE);
    for (int row = 0; row < int(y_size); row++)
    {
      for(int col = 0; col < int(x_size); col++)
      {
        Byte byte = matrix.coeff(row, col) == Scalar(0) ? 255 : 0;
        image_data.GetByte(row, col, 0) = byte;
      }
    }
    hs::imgio::whole::ImageIO image_io;
    image_io.SaveImage(image_path, image_data);
    return 0;
  }
#endif

  Err CompareNormalEquation(const NormalMatrix& normal_matrix,
                            const Gradient& gradient,
                            const MatrixXX& normal_matrix_dense,
                            const VectorX& gradient_dense) const
  {
    Index x_size = normal_matrix_dense.cols();
    Index y_size = normal_matrix_dense.rows();
    const Scalar threshold = Scalar(1e-8);

    Err result = 0;
#ifdef IMAGE_TEST
    MatrixXX normal_diff = MatrixXX::Zero(y_size, x_size);
#endif
    for (Index i = 0; i < y_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        Scalar value = normal_matrix.coeff(i, j);
        Scalar value_dense = normal_matrix_dense(i, j);
        Scalar error_abs = std::abs(value - value_dense);
        Scalar error_rel = error_abs;
        if (value_dense != 0) error_rel /= std::abs(value_dense);
        if (error_abs > threshold && error_rel > threshold)
        {
#ifdef IMAGE_TEST
          normal_diff(i, j) = std::max(error_abs, error_rel);
#endif
          std::cout<<"error_abs:"<<error_abs<<"\n";
          std::cout<<"error_rel:"<<error_rel<<"\n";
          std::cout<<"      normal_matrix["<<i<<","<<j<<"]="<<value<<"\n";
          std::cout<<"normal_matrix_dense["<<i<<","<<j<<"]="<<value_dense<<"\n";
          result = -1;
        }
      }
    }
#ifdef IMAGE_TEST
    SaveMatrixImage(normal_diff, "normal_diff.png");
#endif

    for (Index i = 0; i < x_size; i++)
    {
      Scalar value = gradient[i];
      Scalar value_dense = gradient_dense[i];
      Scalar error_abs = std::abs(value - value_dense);
        Scalar error_rel = error_abs;
        if (value_dense != 0) error_rel /= std::abs(value_dense);
        if (error_abs > threshold && error_rel > threshold)
        {
          std::cout<<"error_abs:"<<error_abs<<"\n";
          std::cout<<"error_rel:"<<error_rel<<"\n";
          std::cout<<"      gradient["<<i<<"]="<<value<<"\n";
          std::cout<<"gradient_dense["<<i<<"]="<<value_dense<<"\n";
          result = -1;
        }
    }

    return result;
  }

private:
  Scalar point_planar_stddev_;
  Scalar point_height_stddev_;
  Scalar image_rotation_stddev_;
  Scalar image_position_stddev_;
  Scalar camera_skew_stddev_;
  Scalar camera_pixel_ratio_stddev_;
};

TEST(TestCameraSharedNormalEquationBuilder, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestCameraSharedNormalEquationBuilder<Scalar>
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
                               intrinsic_params_set);

  VectorFunction vector_function;
  XVector x;
  YVector y_synthetic, y;

  ASSERT_EQ(0, data_generator(vector_function, x, y_synthetic));
  ASSERT_EQ(0, vector_function(x, y));

  Scalar point_planar_stddev = 0.05;
  Scalar point_height_stddev = 0.1;
  Scalar image_rotation_stddev = 0.5;
  Scalar image_position_stddev = 1;
  Scalar camera_skew_stddev = 0.001;
  Scalar camera_pixel_ratio_stddev = 0.001;
  Tester tester(point_planar_stddev,
                point_height_stddev,
                image_rotation_stddev,
                image_position_stddev,
                camera_skew_stddev,
                camera_pixel_ratio_stddev);
  KeyCovariance feature_covariance = KeyCovariance::Identity();
  feature_covariance *= Scalar(4);
  ASSERT_EQ(0, tester.Test(vector_function, x, feature_covariance));
}

}
