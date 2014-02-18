#ifndef _HS_SFM_UNIT_TEST_TEST_CAMERA_SHARED_NORMAL_MLE_META_HPP_
#define _HS_SFM_UNIT_TEST_TEST_CAMERA_SHARED_NORMAL_MLE_META_HPP_

#include "hs_test_utility/test_monte_carlo/normal_mle_meta_forward_declare.hpp"
#include "hs_test_utility/test_monte_carlo/normal_mle_meta_eigen_dense_analytical_x_covariance_calculator_type.hpp"
#include "hs_test_utility/test_monte_carlo/normal_mle_meta_eigen_dense_residuals_calculator_type.hpp"
#include "hs_test_utility/test_monte_carlo/normal_mle_meta_eigen_dense_statistical_x_covariance_calculator_type.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_noised_y_generator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_mahalanobis_distance_calculator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_analytical_x_covariance_calculator.hpp"

namespace hs
{
namespace test
{
namespace mc
{

template <typename _Scalar>
struct NormalMLEXVectorOptimizorType<hs::sfm::ba::CameraSharedVectorFunction<
                                       _Scalar> >
{
  typedef hs::sfm::ba::CameraSharedLevenbergMarquardtOptimizor<
            hs::sfm::ba::CameraSharedVectorFunction<_Scalar> > type;
};

template <typename _Scalar>
struct NormalMLENoisedYVectorGeneratorType<
         EIGEN_VECTOR(_Scalar, Eigen::Dynamic),
         hs::sfm::ba::CameraSharedYCovarianceInverse<_Scalar> >
{
  typedef hs::sfm::ba::CameraSharedNoisedYGenerator<_Scalar> type;
};

template <typename _Scalar>
struct NormalMLEMahalanobisDistanceCalculatorType<
         EIGEN_VECTOR(_Scalar, Eigen::Dynamic)>
{
  typedef hs::sfm::ba::CameraSharedMahalanobisDistanceCalculator<_Scalar> type;
};

template <typename _Scalar>
struct NormalMLEAnalyticalXCovarianceCalculatorType<
         hs::sfm::ba::CameraSharedVectorFunction<_Scalar>,
         EIGEN_VECTOR(_Scalar, Eigen::Dynamic),
         hs::sfm::ba::CameraSharedYCovarianceInverse<_Scalar> >
{
  typedef hs::sfm::ba::CameraSharedAnalyticalXCovarianceCalculator<
    hs::sfm::ba::CameraSharedVectorFunction<_Scalar> > type;
};

}
}
}

#endif
