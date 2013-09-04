#ifndef _UNIT_TEST_TEST_MULTIPLE_VIEW_NORMAL_MLE_META_HPP_
#define _UNIT_TEST_TEST_MULTIPLE_VIEW_NORMAL_MLE_META_HPP_

#include "hs_test_utility/test_monte_carlo/normal_mle_meta_fwd_decl.hpp"
#include "hs_test_utility/test_monte_carlo/normal_mle_meta_eigen_dense.hpp"

#include "hs_sfm/triangulate/multiple_view_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/triangulate/multiple_view_x_covariance_calculator.hpp"
#include "hs_sfm/triangulate/multiple_view_residuals.hpp"

namespace hs
{
namespace test
{
namespace mc
{

template <typename _Scalar>
struct NormalMLEXVectorOptimizorType<
         hs::sfm::triangulate::MultipleViewVectorFunction<_Scalar> >
{
  typedef hs::sfm::triangulate::MultipleViewLevenbergMarquardtOptimizor<
            hs::sfm::triangulate::MultipleViewVectorFunction<_Scalar> > type;
};

template <typename _Scalar>
struct NormalMLEAnalyticXCovarianceCalculatorType<
         hs::sfm::triangulate::MultipleViewVectorFunction<_Scalar> >
{
  typedef hs::sfm::triangulate::MultipleViewXCovarianceCalculator<
            hs::sfm::triangulate::MultipleViewVectorFunction<_Scalar> > type;
};

template <typename _Scalar>
struct NormalMLEResidualsCalculatorType<
         hs::sfm::triangulate::MultipleViewVectorFunction<_Scalar> >
{
  typedef hs::sfm::triangulate::MultipleViewResidualsCalculator<
            hs::sfm::triangulate::MultipleViewVectorFunction<_Scalar> > type;
};

}
}
}

#endif