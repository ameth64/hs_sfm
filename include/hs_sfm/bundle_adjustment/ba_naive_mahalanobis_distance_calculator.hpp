#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_MAHALANOBIS_DISTANCE_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_MAHALANOBIS_DISTANCE_CALCULATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_feature_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BANaiveMahalanobisDistanceCalculator
{
public:
  typedef _Scalar Scalar;
  typedef BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::YVector Residuals;
  typedef BANaiveFeatureCovarianceInverse<Scalar> YCovarianceInverse;

private:
  typedef typename VectorFunction::Index Index;
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;
  typedef typename YCovarianceInverse::Block YCovarianceInverseBlock;

public:
  Scalar operator() (const Residuals& residuals,
                     const YCovarianceInverse& y_covariance_inverse) const
  {
    Index number_of_features = residuals.rows() / 
                               VectorFunction::params_per_feature_;
    if (Index(y_covariance_inverse.blocks.size() != number_of_features))
    {
      return Scalar(0);
    }

    Scalar mahalanobis_distance = Scalar(0);
    for (Index i = 0; i < number_of_features; i++)
    {
      Vector2 feature = residuals.segment(
        i * VectorFunction::params_per_feature_,
        VectorFunction::params_per_feature_);
      const YCovarianceInverseBlock& y_covariance_inverse_block =
        y_covariance_inverse.blocks[i];
      mahalanobis_distance +=
        feature.transpose() * y_covariance_inverse_block *
        feature;
    }

    return mahalanobis_distance;
  }
};

}
}
}

#endif
