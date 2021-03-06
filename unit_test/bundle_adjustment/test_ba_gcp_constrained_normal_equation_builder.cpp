#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/linear_algebra/latraits/vector_eigen.hpp"
#include "hs_math/linear_algebra/lafunc/arithmetic_eigen.hpp"
#include "hs_math/fdjac/forward_finite_difference_sparse_jacobian_matrix_calculator.hpp"

#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_forward_finite_difference_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_noised_y_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_synthetic_data_generator.hpp"

#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_normal_equation_builder.hpp"

namespace
{

template <typename _Scalar>
class TestBAGCPConstrainedNormalEquationBuilder
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef hs::sfm::ba::BAGCPConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;
  typedef hs::sfm::ba::BAGCPConstrainedNormalEquationBuilder<Scalar>
          NormalEquationBuilder;
  typedef typename NormalEquationBuilder::YCovarianceInverse YCovarianceInverse;
  typedef typename NormalEquationBuilder::NormalMatrix NormalMatrix;
  typedef typename NormalEquationBuilder::Gradient Gradient;
  typedef typename NormalEquationBuilder::Residuals Residuals;

  typedef
    hs::math::fdjac::ForwardFiniteDifferenceSparseJacobianMatrixCalculator<
      VectorFunction> FFDJacobianMatrixCalculator;
  typedef typename FFDJacobianMatrixCalculator::JacobianMatrix
                   FFDJacobianMatrix;

  typedef
    hs::sfm::ba::BAGCPConstrainedForwardFiniteDifferenceJacobianMatrixCalculator<
      VectorFunction> BAFFDJacobianMatrixCalculator;
  typedef typename BAFFDJacobianMatrixCalculator::JacobianMatrix
                   BAFFDJacobianMatrix;
  typedef
    hs::sfm::ba::BAGCPConstrainedAnalyticalJacobianMatrixCalculator<
      VectorFunction> BAAnalyticJacobianMatrixCalculator;
  typedef typename BAAnalyticJacobianMatrixCalculator::JacobianMatrix
                   BAAnalyticJacobianMatrix;

  typedef hs::sfm::ba::BAGCPConstrainedNoisedYGenerator<Scalar> NoisedYGenerator;

private:
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic) MatrixXX;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) VectorX;
  typedef EIGEN_MATRIX(Scalar, 2, 2) Matrix22;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

public:
  static Err Test(const VectorFunction& vector_function,
                  const XVector& x,
                  const Matrix22& feature_covariance,
                  const Matrix33& gcp_covariance)
  {
    Err result = 0;

    YVector y;
    if (vector_function(x, y) != 0)
    {
      std::cout<<"Fail to call vector function.\n";
      return -1;
    }

    Index number_of_features = vector_function.number_of_features();
    Index number_of_gcps = vector_function.number_of_gcps();
    MatrixXX dense_y_covariance_inverse;
    YCovarianceInverse y_covariance_inverse;
    if (GenerateYCovarianceInverse(
          feature_covariance,
          gcp_covariance,
          number_of_features,
          number_of_gcps,
          dense_y_covariance_inverse,
          y_covariance_inverse) != 0)
    {
      std::cout<<"GenerateYCovaianceInverse Failed.\n";
      return -1;
    }

    YVector noised_y;
    NoisedYGenerator noised_y_generator;
    if (noised_y_generator(y, y_covariance_inverse, noised_y) != 0)
    {
      std::cout<<"noised y generator failed.\n";
      return -1;
    }

    FFDJacobianMatrixCalculator ffd_jacobian_matrix_calculator(
      Scalar(1e-6), Scalar(1e-10), Scalar(1e-12));
    BAAnalyticJacobianMatrixCalculator ba_analytical_jacobian_matrix_calculator;
    BAFFDJacobianMatrixCalculator ba_ffd_jacobian_matrix_calculator(
      Scalar(1e-6), Scalar(1e-10));

    FFDJacobianMatrix ffd_jacobian_matrix;
    BAAnalyticJacobianMatrix ba_analytical_jacobian_matrix;
    BAFFDJacobianMatrix ba_ffd_jacobian_matrix;

    if (ffd_jacobian_matrix_calculator(vector_function, x,
                                       ffd_jacobian_matrix) != 0)
    {
      std::cout<<"ffd_jacobian_matrix_calculator failed.\n";
      return -1;
    }

    if (ba_analytical_jacobian_matrix_calculator(
      vector_function, x, ba_analytical_jacobian_matrix) != 0)
    {
      std::cout<<"ba_analytical_jacobian_matrix_calculator failed.\n";
      return -1;
    }

    if (ba_ffd_jacobian_matrix_calculator(
      vector_function, x, ba_ffd_jacobian_matrix) != 0)
    {
      std::cout<<"ba_ffd_jacobian_matrix_calculator failed.\n";
      return -1;
    }

    Residuals residuals = y - noised_y;

    NormalMatrix ba_ffd_normal_matrix;
    Gradient ba_ffd_gradient;
    NormalEquationBuilder builder;
    if (builder(ba_ffd_jacobian_matrix,
                residuals,
                y_covariance_inverse,
                ba_ffd_normal_matrix,
                ba_ffd_gradient) != 0)
    {
      std::cout<<"builder failed.\n";
      return -1;
    }

    MatrixXX dense_normal_matrix;
    Index x_size = vector_function.GetXSize();
    dense_normal_matrix.resize(x_size, x_size);
    dense_normal_matrix = ffd_jacobian_matrix.transpose() *
                          dense_y_covariance_inverse *
                          ffd_jacobian_matrix;
    VectorX dense_gradient(x_size);
    dense_gradient = ffd_jacobian_matrix.transpose() *
                     dense_y_covariance_inverse *
                     residuals;

    const Scalar threshold = Scalar(1e-8);
    for (Index i = 0; i < x_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        Scalar dense_value = dense_normal_matrix.coeff(i, j);
        Scalar ba_ffd_value = ba_ffd_normal_matrix.coeff(i, j);
        Scalar err = dense_value - ba_ffd_value;
        if (dense_value != Scalar(0))
        {
          err /= dense_value;
        }
        err = std::abs(err);
        if (err > threshold)
        {
          std::cout<<"dense_normal_matrix["<<i<<", "<<j<<"]:"
                   <<dense_value<<" .\n";
          std::cout<<"ba_ffd_normal_matrix["<<i<<", "<<j<<"]:"
                   <<ba_ffd_value<<" .\n";
          result = -1;
        }
      }
    }

    for (Index i = 0; i < x_size; i++)
    {
      Scalar dense_value = dense_gradient[i];
      Scalar ba_ffd_value = ba_ffd_gradient[i];
      Scalar err = dense_value - ba_ffd_value;
      if (dense_value != Scalar(0))
      {
        err /= dense_value;
      }
      err = std::abs(err);
      if (err > threshold)
      {
        std::cout<<"dense_gradient["<<i<<"]:"
                 <<dense_value<<" .\n";
        std::cout<<"ba_ffd_gradient["<<i<<"]:"
                 <<ba_ffd_value<<" .\n";
        result = -1;
      }
    }

    return result;
  }

private:
  static Err GenerateYCovarianceInverse(
    const Matrix22& feature_covariance,
    const Matrix33& gcp_covariance,
    Index number_of_features,
    Index number_of_gcps,
    MatrixXX& dense_y_covariance_inverse,
    YCovarianceInverse& y_covariance_inverse)
  {
    y_covariance_inverse.naive_y_covariance_inverse.blocks.clear();
    Index feature_size =
      number_of_features * VectorFunction::params_per_feature_;
    Index y_size = feature_size +
                   number_of_gcps * VectorFunction::params_per_point_;
    dense_y_covariance_inverse.resize(y_size, y_size);
    dense_y_covariance_inverse.setZero();
    for (Index i = 0; i < number_of_features; i++)
    {
      dense_y_covariance_inverse.template block<
        VectorFunction::params_per_feature_,
        VectorFunction::params_per_feature_>(
        i * VectorFunction::params_per_feature_,
        i * VectorFunction::params_per_feature_) = feature_covariance.inverse();
      y_covariance_inverse.naive_y_covariance_inverse.blocks.push_back(
        feature_covariance.inverse());
    }
    for (Index i = 0; i < number_of_gcps; i++)
    {
      dense_y_covariance_inverse.template block<
        VectorFunction::params_per_point_,
        VectorFunction::params_per_point_>(
        feature_size + i * VectorFunction::params_per_point_,
        feature_size + i * VectorFunction::params_per_point_) = 
          gcp_covariance.inverse();
      y_covariance_inverse.gcp_blocks.push_back(gcp_covariance.inverse());
    }

    return 0;
  }
};

TEST(TestBAGCPConstrainedNormalEquationBuilder, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef
    hs::sfm::ba::BAGCPConstrainedSyntheticDataGenerator<Scalar, ImageDimension>
    DataGenerator;
  typedef TestBAGCPConstrainedNormalEquationBuilder<Scalar> Test;
  typedef EIGEN_MATRIX(Scalar, 2, 2) Matrix22;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef Test::VectorFunction VectorFunction;
  typedef VectorFunction::XVector XVector;
  typedef VectorFunction::YVector YVector;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 3;
  size_t number_of_cameras_in_strip = 3;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 30;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 1;
  Scalar north_west_angle = 60;
  size_t number_of_gcps = 5;

  VectorFunction vector_function;
  XVector x;
  YVector y;

  DataGenerator data_generator(focal_length_in_metre,
                               number_of_strips,
                               number_of_cameras_in_strip,
                               ground_resolution,
                               image_width,
                               image_height,
                               pixel_size,
                               number_of_points,
                               lateral_overlap_ratio,
                               longitudinal_overlap_ratio,
                               scene_max_height,
                               camera_height_stddev,
                               camera_planar_stddev,
                               camera_rotation_stddev,
                               north_west_angle,
                               number_of_gcps);

  ASSERT_EQ(0, data_generator(vector_function, x, y));
  Scalar focal_length_in_pixel = data_generator.GetFocalLengthInPixel();

  Matrix22 feature_covariance = Matrix22::Identity();
  feature_covariance /= focal_length_in_pixel * focal_length_in_pixel;

  Matrix33 gcp_covariance = Matrix33::Identity();
  gcp_covariance *= Scalar(2);

  ASSERT_EQ(0, Test::Test(vector_function, x,
                          feature_covariance, gcp_covariance));
}

}