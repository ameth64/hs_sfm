#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_MAHALANOBIS_DISTANCE_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_MAHALANOBIS_DISTANCE_CALCULATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_y_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedMahalanobisDistanceCalculator
{
public:
  typedef _Scalar Scalar;
  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::YVector Residuals;
  static const Index params_per_key_ = VectorFunction::params_per_key_;
  typedef CameraSharedYCovarianceInverse<Scalar> YCovarianceInverse;

public:
  Scalar operator() (const Residuals& residuals,
                     const YCovarianceInverse& y_covariance_inverse) const
  {
    size_t number_of_keys = y_covariance_inverse.NumberOfKeys();
    Index key_params_size = y_covariance_inverse.GetKeyParamsSize();
    Index constraints_size = y_covariance_inverse.GetConstraintsSize();

    Scalar mahalanobis_distance = Scalar(0);
    for (size_t i = 0; i < number_of_keys; i++)
    {
      Index key_id = i;
      Scalar key_error =
        residuals.segment(key_id * params_per_key_,
                          params_per_key_).transpose() *
        y_covariance_inverse.GetKeyBlock(i) *
        residuals.segment(key_id * params_per_key_,
                          params_per_key_);
      mahalanobis_distance += key_error;
    }

    for (Index i = 0; i < constraints_size; i++)
    {
      Scalar constraint =
        residuals[key_params_size + i] *
        y_covariance_inverse.GetConstraint(i) *
        residuals[key_params_size + i];
      mahalanobis_distance += constraint;
    }

    return mahalanobis_distance;
  }
};

}
}
}

#endif
