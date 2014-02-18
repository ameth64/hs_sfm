#ifndef _UNIT_TEST_TEST_BA_NAIVE_NORMAL_MLE_META_HPP_
#define _UNIT_TEST_TEST_BA_NAIVE_NORMAL_MLE_META_HPP_

#include "hs_test_utility/test_monte_carlo/normal_mle_meta_forward_declare.hpp"
#include "hs_test_utility/test_monte_carlo/normal_mle_meta_eigen_dense_analytical_x_covariance_calculator_type.hpp"
#include "hs_test_utility/test_monte_carlo/normal_mle_meta_eigen_dense_residuals_calculator_type.hpp"
#include "hs_test_utility/test_monte_carlo/normal_mle_meta_eigen_dense_statistical_x_covariance_calculator_type.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_noised_y_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_mahalanobis_distance_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_analytical_x_covariance_calculator.hpp"

namespace hs
{
namespace test
{
namespace mc
{

template <typename _Scalar>
struct NormalMLEXVectorOptimizorType<hs::sfm::ba::BANaiveVectorFunction<_Scalar> >
{
  typedef hs::sfm::ba::BANaiveLevenbergMarquardtOptimizor<
    hs::sfm::ba::BANaiveVectorFunction<_Scalar> > type;
};

template <typename _Scalar>
struct NormalMLENoisedYVectorGeneratorType<
         EIGEN_VECTOR(_Scalar, Eigen::Dynamic),
         hs::sfm::ba::BANaiveFeatureCovarianceInverse<_Scalar> >
{
  typedef hs::sfm::ba::BANaiveNoisedYGenerator<_Scalar> type;
};

template <typename _Scalar>
struct NormalMLEMahalanobisDistanceCalculatorType<
         EIGEN_VECTOR(_Scalar, Eigen::Dynamic)>
{
  typedef hs::sfm::ba::BANaiveMahalanobisDistanceCalculator<_Scalar> type;
};

template <typename _Scalar>
struct NormalMLEAnalyticalXCovarianceCalculatorType<
  hs::sfm::ba::BANaiveVectorFunction<_Scalar>,
  EIGEN_VECTOR(_Scalar, Eigen::Dynamic),
  hs::sfm::ba::BANaiveFeatureCovarianceInverse<_Scalar> >
{
  typedef hs::sfm::ba::BANaiveAnalyticalXCovarianceCalculator<
    hs::sfm::ba::BANaiveVectorFunction<_Scalar> > type;
};

}
}
}

#endif
