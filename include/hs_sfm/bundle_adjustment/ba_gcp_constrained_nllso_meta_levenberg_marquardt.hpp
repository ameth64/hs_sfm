#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_NLLSO_LEVENBERG_MARQUARDT_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_NLLSO_LEVENBERG_MARQUARDT_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_optimizor/nllso/meta_forward_declare.hpp"
#include "hs_optimizor/nllso/meta_eigen_residuals_calculator_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_x_vector_updater_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_x_vector_norm_calculator_type.hpp"
#include "hs_optimizor/nllso/meta_eigen_delta_x_vector_norm_calculator_type.hpp"

#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_y_covariance_inverse.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_jacobian_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_gradient.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_forward_finite_difference_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_normal_equation_builder.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_augmentor.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_matrix_max_diagonal_value_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_equation_solver.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_gradient_delta_x_dot_product_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_gradient_infinite_norm_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_mahalanobis_distance_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_nllso_meta_levenberg_marquardt.hpp"

namespace hs
{
namespace optimizor
{
namespace nllso
{

#ifdef BA_GCP_CONSTRAINED_VECTOR_FUNCTION
#undef BA_GCP_CONSTRAINED_VECTOR_FUNCTION
#endif
#define BA_GCP_CONSTRAINED_VECTOR_FUNCTION \
  hs::sfm::ba::BAGCPConstrainedVectorFunction<_Scalar>

#ifdef BA_GCP_CONSTRAINED_JACOBIAN_MATRIX
#undef BA_GCP_CONSTRAINED_JACOBIAN_MATRIX
#endif
#define BA_GCP_CONSTRAINED_JACOBIAN_MATRIX \
  hs::sfm::ba::BAGCPConstrainedJacobianMatrix<\
    _Scalar,\
    typename BA_GCP_CONSTRAINED_VECTOR_FUNCTION::Index,\
    BA_GCP_CONSTRAINED_VECTOR_FUNCTION::params_per_feature_,\
    BA_GCP_CONSTRAINED_VECTOR_FUNCTION::params_per_camera_,\
    BA_GCP_CONSTRAINED_VECTOR_FUNCTION::params_per_point_>

template <typename _Scalar>
struct JacobianMatrixCalculatorType<
  BA_GCP_CONSTRAINED_VECTOR_FUNCTION >
{
  typedef BA_GCP_CONSTRAINED_VECTOR_FUNCTION VectorFunction;
  typedef hs::sfm::ba::BAGCPConstrainedAnalyticalJacobianMatrixCalculator<
            VectorFunction> type;
};

template <typename _Scalar>
struct YCovarianceInverseType<
  BA_GCP_CONSTRAINED_VECTOR_FUNCTION >
{
  typedef hs::sfm::ba::BAGCPConstrainedYCovarianceInverse<_Scalar> type;
};

template <typename _Scalar>
struct NormalEquationBuilderType<
  EIGEN_VECTOR(_Scalar, Eigen::Dynamic),
  BA_GCP_CONSTRAINED_JACOBIAN_MATRIX,
  hs::sfm::ba::BAGCPConstrainedYCovarianceInverse<_Scalar> >
{
  typedef hs::sfm::ba::BAGCPConstrainedNormalEquationBuilder<_Scalar> type;
};

template <typename _Scalar>
struct MahalanobisDistanceCalculatorType<
  hs::sfm::ba::BAGCPConstrainedYCovarianceInverse<_Scalar>,
  EIGEN_VECTOR(_Scalar, Eigen::Dynamic) >
{
  typedef hs::sfm::ba::BAGCPConstrainedMahalanobisDistanceCalculator<_Scalar>
          type;
};

}
}
}

#endif