#ifndef _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_MAHALANOBIS_DISTANCE_CALCULATOR_HPP_
#define _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_MAHALANOBIS_DISTANCE_CALCULATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/projective/mle/intrinsic_constrained_vector_function.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_y_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

template <typename _Scalar>
class IntrinsicConstrainedMahalanobisDistanceCalculator
{
public:
  typedef _Scalar Scalar;
  typedef IntrinsicConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::YVector Residuals;
  typedef IntrinsicConstrainedYCovarianceInverse<Scalar> YCovarianceInverse;

private:
  typedef typename VectorFunction::Index Index;
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;
  typedef typename YCovarianceInverse::KeyBlock KeyBlock;

public:
  Scalar operator() (const Residuals& residuals,
                     const YCovarianceInverse& y_covariance_inverse) const
  {
    Index number_of_keys = Index(y_covariance_inverse.key_blocks.size());

    Scalar mahalanobis_distance = Scalar(0);
    for (Index i = 0; i < number_of_keys; i++)
    {
      Vector2 key = residuals.segment(i * 2, 2);
      const KeyBlock& key_block = y_covariance_inverse.key_blocks[i];
      mahalanobis_distance += key.transpose() * key_block * key;
    }
    mahalanobis_distance +=
      residuals[number_of_keys * 2 + 0] * residuals[number_of_keys * 2 + 0] /
      y_covariance_inverse.focal_length_stddev /
      y_covariance_inverse.focal_length_stddev +
      residuals[number_of_keys * 2 + 1] * residuals[number_of_keys * 2 + 1] /
      y_covariance_inverse.skew_stddev /
      y_covariance_inverse.skew_stddev +
      residuals[number_of_keys * 2 + 2] * residuals[number_of_keys * 2 + 2] /
      y_covariance_inverse.principal_point_x_stddev /
      y_covariance_inverse.principal_point_x_stddev +
      residuals[number_of_keys * 2 + 3] * residuals[number_of_keys * 2 + 3] /
      y_covariance_inverse.principal_point_y_stddev /
      y_covariance_inverse.principal_point_y_stddev +
      residuals[number_of_keys * 2 + 4] * residuals[number_of_keys * 2 + 4] /
      y_covariance_inverse.pixel_ratio_stddev /
      y_covariance_inverse.pixel_ratio_stddev;

      return mahalanobis_distance;
  }
};

}
}
}

#endif