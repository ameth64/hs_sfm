#ifndef _UNIT_TEST_TEST_BA_NAIVE_NORMAL_MLE_META_HPP_
#define _UNIT_TEST_TEST_BA_NAIVE_NORMAL_MLE_META_HPP_

#include "hs_test_utility/test_monte_carlo/normal_mle_meta_fwd_decl.hpp"
#include "hs_test_utility/test_monte_carlo/normal_mle_meta_eigen_dense.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_xcov_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_residuals.hpp"

#include "test_ba_naive_noised_y_generator.hpp"

namespace hs
{
namespace test
{
namespace mc
{

template <typename _Scalar>
struct NormalMLEXVectorOptimizorType<hs::sfm::ba::BANaiveVecFunc<_Scalar> >
{
  typedef hs::sfm::ba::BANaiveLevenbergMarquardtOptimizor<
    hs::sfm::ba::BANaiveVecFunc<_Scalar> > type;
};

template <typename _Scalar>
struct NormalMLEAnalyticXCovarianceCalculatorType<
         hs::sfm::ba::BANaiveVecFunc<_Scalar> >
{
  typedef hs::sfm::ba::BANaiveXCovCalculator<
    hs::sfm::ba::BANaiveVecFunc<_Scalar> > type;
};

template <typename _Scalar>
struct NormalMLENoisedYVectorGeneratorType<
         EIGEN_VEC(_Scalar, Eigen::Dynamic),
         hs::sfm::ba::BANaiveFeatCovInv<_Scalar> >
{
  typedef hs::sfm::ba::BANaiveNoisedYGenerator<_Scalar> type;
};

template <typename _Scalar>
struct NormalMLEResidualsCalculatorType<
         hs::sfm::ba::BANaiveVecFunc<_Scalar> >
{
  typedef hs::sfm::ba::BANaiveResidualsCalc<
    hs::sfm::ba::BANaiveVecFunc<_Scalar> > type;
};

}
}
}

#endif