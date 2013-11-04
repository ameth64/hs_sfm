#include <iostream>

#include <gtest/gtest.h>

#include "hs_sfm/projective/mle/intrinsic_constrained_synthetic_generator.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_noised_y_generator.hpp"

#include "hs_sfm/projective/mle/intrinsic_constrained_normal_equation_builder.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestIntrinsicConstrainedNormalEquationBuilder
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
  typedef hs::sfm::projective::IntrinsicConstrainedSyntheticGenerator<
            Scalar, ImageDimension> SyntheticGenerator;
  typedef hs::sfm::projective::IntrinsicConstrainedVectorFunction<Scalar>
          VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef
    hs::sfm::projective::IntrinsicConstrainedAnalyticalJacobianMatrixCalculator<
      VectorFunction> JacobianMatrixCalculator;
  typedef typename JacobianMatrixCalculator::JacobianMatrix
          JacobianMatrix;
  typedef hs::sfm::projective::IntrinsicConstrainedNormalEquationBuilder<Scalar>
          NormalEquationBuilder;
  typedef typename NormalEquationBuilder::Residuals Residuals;
  typedef typename NormalEquationBuilder::NormalMatrix NormalMatrix;
  typedef typename NormalEquationBuilder::Gradient Gradient;
  typedef typename NormalEquationBuilder::YCovarianceInverse YCovarianceInverse;

  typedef hs::sfm::projective::IntrinsicConstrainedNoisedYGenerator<Scalar>
          NoisedYGenerator;
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic)
          DenseYCovarianceInverse;

public:
  TestIntrinsicConstrainedNormalEquationBuilder(
    Scalar focal_length_in_metre,
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    size_t number_of_points,
    Scalar scene_max_height,
    Scalar camera_height_stddev,
    Scalar camera_planar_stddev,
    Scalar camera_rot_stddev,
    Scalar north_west_angle,
    Scalar skew,
    Scalar principal_point_x,
    Scalar principal_point_y,
    Scalar pixel_ratio,
    Scalar key_stddev,
    Scalar focal_length_stddev,
    Scalar skew_stddev,
    Scalar principal_point_x_stddev,
    Scalar principal_point_y_stddev,
    Scalar pixel_ratio_stddev)
    : synthetic_generator_(focal_length_in_metre,
                           ground_resolution,
                           image_width,
                           image_height,
                           pixel_size,
                           number_of_points,
                           scene_max_height,
                           camera_height_stddev,
                           camera_planar_stddev,
                           camera_rot_stddev,
                           north_west_angle,
                           skew,
                           principal_point_x,
                           principal_point_y,
                           pixel_ratio),
      key_stddev_(key_stddev),
      focal_length_stddev_(focal_length_stddev),
      skew_stddev_(skew_stddev),
      principal_point_x_stddev_(principal_point_x_stddev),
      principal_point_y_stddev_(principal_point_y_stddev),
      pixel_ratio_stddev_(pixel_ratio_stddev) {}

  Err Test()
  {
    VectorFunction vector_function;
    XVector x;
    YVector y;
    if (synthetic_generator_(vector_function, x, y) != 0)
    {
      return -1;
    }

    JacobianMatrixCalculator jacobian_matrix_calculator;
    JacobianMatrix jacobian_matrix;
    if (jacobian_matrix_calculator(vector_function, x, jacobian_matrix) != 0)
    {
      return -1;
    }

    Index number_of_points = Index(vector_function.points().size());
    YCovarianceInverse y_covariance_inverse;
    DenseYCovarianceInverse dense_y_covariance_inverse;
    GenerateYCovarianceInverse(number_of_points,
                               y_covariance_inverse,
                               dense_y_covariance_inverse);

    NoisedYGenerator noised_y_generator;
    YVector noised_y;
    if (noised_y_generator(y, y_covariance_inverse, noised_y) != 0)
    {
      return -1;
    }

    Residuals residuals = y - noised_y;

    NormalEquationBuilder builder;
    NormalMatrix normal_matrix;
    Gradient gradient;
    if (builder(jacobian_matrix, residuals, y_covariance_inverse,
                normal_matrix, gradient) != 0)
    {
      return -1;
    }

    NormalMatrix dense_normal_matrix = jacobian_matrix.transpose() *
                                       dense_y_covariance_inverse *
                                       jacobian_matrix;
    Gradient dense_gradient = jacobian_matrix.transpose() *
                              dense_y_covariance_inverse *
                              residuals;

    Index x_size = vector_function.GetXSize();
    const Scalar threshold = Scalar(1e-8);
    Err result = 0;
    for (Index i = 0; i < x_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        Scalar dense_value = dense_normal_matrix.coeff(i, j);
        Scalar analytical_value = normal_matrix.coeff(i, j);
        Scalar err = dense_value - analytical_value;
        if (dense_value != Scalar(0))
        {
          err /= dense_value;
        }
        err = std::abs(err);
        if (err > threshold)
        {
          std::cout<<"dense_normal_matrix["<<i<<", "<<j<<"]:"
                   <<dense_value<<" .\n";
          std::cout<<"normal_matrix["<<i<<", "<<j<<"]:"
                   <<analytical_value<<" .\n";
          result = -1;
        }
      }
    }

    return result;
  }

private:
  void GenerateYCovarianceInverse(
    Index number_of_points,
    YCovarianceInverse& y_covariance_inverse,
    DenseYCovarianceInverse& dense_y_covariance_inverse) const
  {
    typedef typename YCovarianceInverse::KeyBlock KeyBlock;
    Index y_size = number_of_points * 2 + 5;
    dense_y_covariance_inverse.resize(y_size, y_size);
    dense_y_covariance_inverse.setZero();
    KeyBlock key_block = KeyBlock::Identity() / key_stddev_ / key_stddev_;
    for (Index i = 0; i < number_of_points; i++)
    {
      y_covariance_inverse.key_blocks.push_back(key_block);
      dense_y_covariance_inverse.template block<2, 2>(i * 2, i * 2) = key_block;
    }

    y_covariance_inverse.focal_length_stddev = focal_length_stddev_;
    y_covariance_inverse.skew_stddev = skew_stddev_;
    y_covariance_inverse.principal_point_x_stddev = principal_point_x_stddev_;
    y_covariance_inverse.principal_point_y_stddev = principal_point_y_stddev_;
    y_covariance_inverse.pixel_ratio_stddev = pixel_ratio_stddev_;

    dense_y_covariance_inverse(number_of_points * 2 + 0,
                               number_of_points * 2 + 0) =
      1 / focal_length_stddev_ / focal_length_stddev_;
    dense_y_covariance_inverse(number_of_points * 2 + 1,
                               number_of_points * 2 + 1) =
      1 / skew_stddev_ / skew_stddev_;
    dense_y_covariance_inverse(number_of_points * 2 + 2,
                               number_of_points * 2 + 2) =
      1 / principal_point_x_stddev_ / principal_point_x_stddev_;
    dense_y_covariance_inverse(number_of_points * 2 + 3,
                               number_of_points * 2 + 3) =
      1 / principal_point_y_stddev_ / principal_point_y_stddev_;
    dense_y_covariance_inverse(number_of_points * 2 + 4,
                               number_of_points * 2 + 4) =
      1 / pixel_ratio_stddev_ / pixel_ratio_stddev_;
  }

private:
  SyntheticGenerator synthetic_generator_;
  Scalar key_stddev_;
  Scalar focal_length_stddev_;
  Scalar skew_stddev_;
  Scalar principal_point_x_stddev_;
  Scalar principal_point_y_stddev_;
  Scalar pixel_ratio_stddev_;
};

TEST(TestIntrinsicConstrainedNormalEquationBuilder, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestIntrinsicConstrainedNormalEquationBuilder<
            Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 10;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rot_stddev = 1;
  Scalar north_west_angle = 60;
  Scalar skew = 0.0001;
  Scalar principal_point_x = 20.0;
  Scalar principal_point_y = 30.0;
  Scalar pixel_ratio = 1.0001;
  Scalar key_stddev = 1.0;
  Scalar focal_length_stddev = 1;
  Scalar skew_stddev = 0.00001;
  Scalar principal_point_x_stddev = 1;
  Scalar principal_point_y_stddev = 1;
  Scalar pixel_ratio_stddev = 0.00001;

  Test test(focal_length_in_metre,
            ground_resolution,
            image_width,
            image_height,
            pixel_size,
            number_of_points,
            scene_max_height,
            camera_height_stddev,
            camera_planar_stddev,
            camera_rot_stddev,
            north_west_angle,
            skew,
            principal_point_x,
            principal_point_y,
            pixel_ratio,
            key_stddev,
            focal_length_stddev,
            skew_stddev,
            principal_point_x_stddev,
            principal_point_y_stddev,
            pixel_ratio_stddev);

  ASSERT_EQ(0, test.Test());
}

}