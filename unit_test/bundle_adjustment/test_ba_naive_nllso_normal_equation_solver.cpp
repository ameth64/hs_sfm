#include <iostream>

#include <gtest/gtest.h>

#include "hs_sfm/bundle_adjustment/ba_naive_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_equation_builder.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_noised_y_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_synthetic_data_generator.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_normal_matrix_max_diagonal_value_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_augmentor.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_equation_solver.hpp"

namespace
{

template <typename _Scalar>
class TestBANaiveNLLSONormalEquationSolver
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef hs::sfm::ba::BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;
  typedef hs::sfm::ba::BANaiveNormalEquationBuilder<Scalar>
          NormalEquationBuilder;
  typedef typename NormalEquationBuilder::YCovarianceInverse
                   YCovarianceInverse;
  typedef typename NormalEquationBuilder::NormalMatrix NormalMatrix;
  typedef typename NormalEquationBuilder::Gradient Gradient;
  typedef YVector Residuals;

  typedef hs::sfm::ba::BANaiveAnalyticalJacobianMatrixCalculator<
            VectorFunction>
          BAAnalyticJacobianMatrixCalculator;
  typedef typename BAAnalyticJacobianMatrixCalculator::JacobianMatrix
          BAAnalyticJacobianMatrix;

  typedef hs::sfm::ba::BANaiveNormalEquationSolver<Scalar>
          NormalEquationSolver;
  typedef typename NormalEquationSolver::DeltaXVector DeltaXVector;

  typedef hs::sfm::ba::BANaiveAugmentor<Scalar> Augmentor;
  typedef typename Augmentor::AugmentedNormalMatrix AugmentedNormalMatrix;

  typedef hs::sfm::ba::BANaiveNoisedYGenerator<Scalar> NoisedYGenerator;

  typedef hs::sfm::ba::BANaiveNormalMatrixMaxDiagonalValueCalculator<Scalar>
          MaxDiagonalValueCalculator;
  
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic)
          MatrixXX;
  typedef EIGEN_MATRIX(Scalar, 2, 2) Matrix22;

  static Err Test(const VectorFunction& vector_function,
                  const XVector& x,
                  const Matrix22& feature_covariance)
  {
    Err result = 0;

    YVector y;
    if (vector_function(x, y) != 0)
    {
      std::cout<<"Fail to call vector function.\n";
      return -1;
    }

    Index number_of_features = vector_function.number_of_features();
    MatrixXX dense_y_covariance_inverse;
    YCovarianceInverse y_covariance_inverse;
    if (GenerateFeatureCovaianceInverse(feature_covariance,
                                        number_of_features,
                                        dense_y_covariance_inverse,
                                        y_covariance_inverse) != 0)
    {
      std::cout<<"GenerateFeatureCovaianceInverse Failed.\n";
      return -1;
    }

    YVector noised_y;
    NoisedYGenerator noised_y_generator;
    if (noised_y_generator(y, y_covariance_inverse, noised_y) != 0)
    {
      std::cout<<"noised_y_generator failed.\n";
      return -1;
    }

    BAAnalyticJacobianMatrixCalculator ba_analytical_jacobian_matrix_calculator;
    BAAnalyticJacobianMatrix ba_analytical_jacobian_matrix;
    if (ba_analytical_jacobian_matrix_calculator(
          vector_function, x, ba_analytical_jacobian_matrix) != 0)
    {
      std::cout<<"ba_analytical_jacobian_matrix_calculator failed.\n";
      return -1;
    }

    Residuals residuals = y - noised_y;

    NormalMatrix normal_matrix;
    Gradient gradient;
    NormalEquationBuilder builder;
    if (builder(ba_analytical_jacobian_matrix,
                residuals, y_covariance_inverse,
                normal_matrix, gradient) != 0)
    {
      std::cout<<"builder failed.\n";
      return -1;
    }

    MaxDiagonalValueCalculator max_diagonal_value_calculator;
    Scalar max_diagonal_value =
      max_diagonal_value_calculator(normal_matrix);

    Scalar mu = /*Scalar(1e-3) * */max_diagonal_value;

    Augmentor augmentor;
    AugmentedNormalMatrix augmented_normal_matrix =
      augmentor(normal_matrix, mu);

    NormalEquationSolver solver;
    DeltaXVector delta_x;
    if (solver(augmented_normal_matrix, gradient, delta_x) != 0)
    {
      std::cout<<"solver failed.\n";
      return -1;
    }

    MatrixXX dense_normal_matrix;
    XVector dense_gradient;
    Index x_size = vector_function.GetXSize();
    dense_normal_matrix.resize(x_size, x_size);
    dense_gradient.resize(x_size);
    for (Index i = 0; i < x_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        dense_normal_matrix(i, j) = normal_matrix.coeff(i, j);
      }
      dense_gradient[i] = gradient[i];
      dense_normal_matrix(i, i) += mu;
    }

    Scalar threshold = Scalar(1e-8);
    
    if ((dense_normal_matrix * delta_x).isApprox(dense_gradient, threshold))
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }

private:
  static Err GenerateFeatureCovaianceInverse(
    const Matrix22& feature_covariance,
    Index number_of_features,
    MatrixXX& dense_y_covariance_inverse,
    YCovarianceInverse& y_covariance_inverse)
  {
    y_covariance_inverse.blocks.clear();
    dense_y_covariance_inverse.resize(
      number_of_features * VectorFunction::params_per_feature_,
      number_of_features * VectorFunction::params_per_feature_);
    dense_y_covariance_inverse.setZero();
    for (Index i = 0; i < number_of_features; i++)
    {
      dense_y_covariance_inverse.template block<
        VectorFunction::params_per_feature_,
        VectorFunction::params_per_feature_>(
        i * VectorFunction::params_per_feature_,
        i * VectorFunction::params_per_feature_) = feature_covariance.inverse();
      y_covariance_inverse.blocks.push_back(feature_covariance.inverse());
    }

    return 0;
  }
};

TEST(TestBANaiveNLLSONormalEquationSolver, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef hs::sfm::ba::BANaiveSyntheticDataGenerator<Scalar, ImageDimension>
          DataGenerator;
  typedef TestBANaiveNLLSONormalEquationSolver<Scalar> Test;
  typedef Test::Matrix22 Matrix22;
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
                               north_west_angle);
  ASSERT_EQ(0, data_generator(vector_function, x, y));
  Scalar focal_length_in_pixel = data_generator.GetFocalLengthInPixel();

  Matrix22 feature_covariance = Matrix22::Identity();
  feature_covariance /= focal_length_in_pixel * focal_length_in_pixel;

  ASSERT_EQ(0, Test::Test(vector_function, x, feature_covariance));
}

}
