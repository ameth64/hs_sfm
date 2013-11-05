#ifndef _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_NLLSO_META_LEVENBERG_MARQUARDT_HPP_
#define _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_NLLSO_META_LEVENBERG_MARQUARDT_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_optimizor/nllso/meta_forward_declare.hpp"
#include "hs_optimizor/nllso/meta_eigen_residuals_calculator_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_x_vector_updater_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_x_vector_norm_calculator_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_delta_x_vector_norm_calculator_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_normal_matrix_augmentor_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_normal_equation_solver_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_normal_matrix_max_diagonal_value_calculator_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_gradient_delta_x_dot_product_calculator_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_gradient_infinite_norm_calculator_type.hpp"

#include "hs_sfm/projective/mle/intrinsic_constrained_vector_function.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_y_covariance_inverse.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_normal_equation_builder.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_mahalanobis_distance_calculator.hpp"

namespace hs
{
namespace optimizor
{
namespace nllso
{

#ifdef INTRINSIC_CONSTRAINED_VECTOR_FUNCTION
#undef INTRINSIC_CONSTRAINED_VECTOR_FUNCTION
#endif
#define INTRINSIC_CONSTRAINED_VECTOR_FUNCTION \
  hs::sfm::projective::IntrinsicConstrainedVectorFunction<_Scalar>

template <typename _Scalar>
struct JacobianMatrixCalculatorType<
  INTRINSIC_CONSTRAINED_VECTOR_FUNCTION>
{
  typedef INTRINSIC_CONSTRAINED_VECTOR_FUNCTION VectorFunction;
  typedef
    hs::sfm::projective::IntrinsicConstrainedAnalyticalJacobianMatrixCalculator<
      VectorFunction> type;
};

template <typename _Scalar>
struct YCovarianceInverseType<
  INTRINSIC_CONSTRAINED_VECTOR_FUNCTION>
{
  typedef hs::sfm::projective::IntrinsicConstrainedYCovarianceInverse<_Scalar>
          type;
};

template <typename _Scalar>
struct NormalEquationBuilderType<
  EIGEN_VECTOR(_Scalar, Eigen::Dynamic),
  EIGEN_MATRIX(_Scalar, Eigen::Dynamic, Eigen::Dynamic),
  hs::sfm::projective::IntrinsicConstrainedYCovarianceInverse<_Scalar> >
{
  typedef
    hs::sfm::projective::IntrinsicConstrainedNormalEquationBuilder<_Scalar>
      type;
};

template <typename _Scalar>
struct MahalanobisDistanceCalculatorType<
  hs::sfm::projective::IntrinsicConstrainedYCovarianceInverse<_Scalar>,
  EIGEN_VECTOR(_Scalar, Eigen::Dynamic)>
{
  typedef
    hs::sfm::projective::IntrinsicConstrainedMahalanobisDistanceCalculator<
      _Scalar> type;
};

}
}
}

#endif