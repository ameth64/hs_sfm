#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NLLSO_META_LEVENBERG_MARQUARDT_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NLLSO_META_LEVENBERG_MARQUARDT_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_optimizor/nllso/meta_forward_declare.hpp"
#include "hs_optimizor/nllso/meta_eigen_residuals_calculator_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_x_vector_updater_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_x_vector_norm_calculator_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_delta_x_vector_norm_calculator_type.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_feature_covariance_inverse.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_jacobian_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_gradient.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_forward_finite_difference_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_equation_builder.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_augmentor.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_matrix_max_diagonal_value_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_equation_solver.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_gradient_delta_x_dot_product_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_gradient_infinite_norm_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_mahalanobis_distance_calculator.hpp"


namespace hs
{
namespace optimizor
{
namespace nllso
{

#ifdef BA_NAIVE_VECTOR_FUNCTION
#undef BA_NAIVE_VECTOR_FUNCTION
#endif
#define BA_NAIVE_VECTOR_FUNCTION hs::sfm::ba::BANaiveVectorFunction<_Scalar>

#ifdef BA_NAIVE_JACOBIAN_MATRIX
#undef BA_NAIVE_JACOBIAN_MATRIX
#endif
#define BA_NAIVE_JACOBIAN_MATRIX \
  hs::sfm::ba::BANaiveJacobianMatrix<\
    _Scalar,\
    typename BA_NAIVE_VECTOR_FUNCTION::Index,\
    BA_NAIVE_VECTOR_FUNCTION::params_per_feature_,\
    BA_NAIVE_VECTOR_FUNCTION::params_per_camera_,\
    BA_NAIVE_VECTOR_FUNCTION::params_per_point_>

#ifdef BA_NAIVE_AUGMENTED_NORMAL_MATRIX
#undef BA_NAIVE_AUGMENTED_NORMAL_MATRIX
#endif
#define BA_NAIVE_AUGMENTED_NORMAL_MATRIX \
  hs::sfm::ba::BANaiveAugmentedNormalMatrix<\
    _Scalar,\
    typename BA_NAIVE_VECTOR_FUNCTION::Index,\
    BA_NAIVE_VECTOR_FUNCTION::params_per_camera_,\
    BA_NAIVE_VECTOR_FUNCTION::params_per_point_>

#ifdef BA_NAIVE_NORMAL_MATRIX
#undef BA_NAIVE_NORMAL_MATRIX
#endif
#define BA_NAIVE_NORMAL_MATRIX \
  hs::sfm::ba::BANaiveNormalMatrix<\
    _Scalar,\
    typename BA_NAIVE_VECTOR_FUNCTION::Index,\
    BA_NAIVE_VECTOR_FUNCTION::params_per_camera_,\
    BA_NAIVE_VECTOR_FUNCTION::params_per_point_>

#ifdef BA_NAIVE_GRADIENT
#undef BA_NAIVE_GRADIENT
#endif
#define BA_NAIVE_GRADIENT \
  hs::sfm::ba::BANaiveGradient<\
    _Scalar,\
    typename BA_NAIVE_VECTOR_FUNCTION::Index,\
    BA_NAIVE_VECTOR_FUNCTION::params_per_camera_,\
    BA_NAIVE_VECTOR_FUNCTION::params_per_point_>

template <typename _Scalar>
struct JacobianMatrixCalculatorType<
  BA_NAIVE_VECTOR_FUNCTION >
{
  typedef BA_NAIVE_VECTOR_FUNCTION VectorFunction;
  typedef hs::sfm::ba::BANaiveAnalyticalJacobianMatrixCalculator<
            VectorFunction> type;
};

template <typename _Scalar>
struct YCovarianceInverseType<
  BA_NAIVE_VECTOR_FUNCTION>
{
  typedef hs::sfm::ba::BANaiveFeatureCovarianceInverse<_Scalar> type;
};

template <typename _Scalar>
struct NormalEquationBuilderType<
  EIGEN_VECTOR(_Scalar, Eigen::Dynamic),
  BA_NAIVE_JACOBIAN_MATRIX,
  hs::sfm::ba::BANaiveFeatureCovarianceInverse<_Scalar> >
{
  typedef hs::sfm::ba::BANaiveNormalEquationBuilder<_Scalar> type;
};

template <typename _Scalar>
struct NormalEquationSolverType<
  BA_NAIVE_AUGMENTED_NORMAL_MATRIX,
  BA_NAIVE_GRADIENT>
{
  typedef hs::sfm::ba::BANaiveNormalEquationSolver<_Scalar> type;
};

template <typename _Scalar>
struct NormalMatrixAugmentorType<
  BA_NAIVE_NORMAL_MATRIX>
{
  typedef hs::sfm::ba::BANaiveAugmentor<_Scalar> type;
};

template <typename _Scalar>
struct NormalMatrixMaxDiagonalValueCalculatorType<
  BA_NAIVE_NORMAL_MATRIX>
{
  typedef hs::sfm::ba::BANaiveNormalMatrixMaxDiagonalValueCalculator<_Scalar>
          type;
};

template <typename _Scalar>
struct GradientDeltaXDotProductCalculatorType<
  BA_NAIVE_GRADIENT,
  EIGEN_VECTOR(_Scalar, Eigen::Dynamic)>
{
  typedef hs::sfm::ba::BANaiveGradientDeltaXDotProductCalculator<_Scalar> type;
};

template <typename _Scalar>
struct MahalanobisDistanceCalculatorType<
  hs::sfm::ba::BANaiveFeatureCovarianceInverse<_Scalar>,
  EIGEN_VECTOR(_Scalar, Eigen::Dynamic)>
{
  typedef hs::sfm::ba::BANaiveMahalanobisDistanceCalculator<_Scalar> type;
};

template <typename _Scalar>
struct GradientInfiniteNormCalculatorType<
  BA_NAIVE_GRADIENT>
{
  typedef hs::sfm::ba::BANaiveGradientInfiniteNormCalculator<_Scalar> type;
};

}
}

}

#endif
