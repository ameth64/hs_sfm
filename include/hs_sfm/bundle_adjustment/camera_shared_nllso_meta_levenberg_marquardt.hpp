#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NLLSO_META_LEVENBERG_MARQUARDT_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NLLSO_META_LEVENBERG_MARQUARDT_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_optimizor/nllso/meta_forward_declare.hpp"
#include "hs_optimizor/nllso/meta_eigen_residuals_calculator_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_x_vector_updater_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_x_vector_norm_calculator_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_delta_x_vector_norm_calculator_type.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_y_covariance_inverse.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_normal_matrix.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_jacobian_matrix.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_gradient.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_forward_finite_difference_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_normal_equation_builder.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_augmentor.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_normal_matrix_max_diagonal_value_calculator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_normal_equation_solver.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_gradient_delta_x_dot_product_calculator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_gradient_infinite_norm_calculator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_mahalanobis_distance_calculator.hpp"

namespace hs
{
namespace optimizor
{
namespace nllso
{

template <typename _Scalar>
struct JacobianMatrixCalculatorType<
  hs::sfm::ba::CameraSharedVectorFunction<_Scalar> >
{
  typedef hs::sfm::ba::CameraSharedVectorFunction<_Scalar> VectorFunction;
  typedef hs::sfm::ba::CameraSharedAnalyticalJacobianMatrixCalculator<
            VectorFunction> type;
};

template <typename _Scalar>
struct YCovarianceInverseType<
  hs::sfm::ba::CameraSharedVectorFunction<_Scalar> >
{
  typedef hs::sfm::ba::CameraSharedYCovarianceInverse<_Scalar> type;
};

template <typename _Scalar>
struct NormalEquationBuilderType<
  EIGEN_VECTOR(_Scalar, Eigen::Dynamic),
  hs::sfm::ba::CameraSharedJacobianMatrix<_Scalar>,
  hs::sfm::ba::CameraSharedYCovarianceInverse<_Scalar> >
{
  typedef hs::sfm::ba::CameraSharedNormalEquationBuilder<_Scalar> type;
};

template <typename _Scalar>
struct NormalEquationSolverType<
  hs::sfm::ba::CameraSharedAugmentedNormalMatrix<_Scalar>,
  hs::sfm::ba::CameraSharedGradient<_Scalar> >
{
  typedef hs::sfm::ba::CameraSharedNormalEquationSolver<_Scalar> type;
};

template <typename _Scalar>
struct NormalMatrixAugmentorType<
  hs::sfm::ba::CameraSharedNormalMatrix<_Scalar> >
{
  typedef hs::sfm::ba::CameraSharedAugmentor<_Scalar> type;
};

template <typename _Scalar>
struct NormalMatrixMaxDiagonalValueCalculatorType<
  hs::sfm::ba::CameraSharedNormalMatrix<_Scalar> >
{
  typedef
    hs::sfm::ba::CameraSharedNormalMatrixMaxDiagonalValueCalculator<_Scalar>
    type;
};

template <typename _Scalar>
struct GradientDeltaXDotProductCalculatorType<
  hs::sfm::ba::CameraSharedGradient<_Scalar>,
  EIGEN_VECTOR(_Scalar, Eigen::Dynamic)>
{
  typedef hs::sfm::ba::CameraSharedGradientDeltaXDotProductCalculator<_Scalar>
          type;
};

template <typename _Scalar>
struct MahalanobisDistanceCalculatorType<
  hs::sfm::ba::CameraSharedYCovarianceInverse<_Scalar>,
  EIGEN_VECTOR(_Scalar, Eigen::Dynamic)>
{
  typedef hs::sfm::ba::CameraSharedMahalanobisDistanceCalculator<_Scalar>
          type;
};

template <typename _Scalar>
struct GradientInfiniteNormCalculatorType<
  hs::sfm::ba::CameraSharedGradient<_Scalar> >
{
  typedef hs::sfm::ba::CameraSharedGradientInfiniteNormCalculator<_Scalar>
          type;
};

}
}
}

#endif
