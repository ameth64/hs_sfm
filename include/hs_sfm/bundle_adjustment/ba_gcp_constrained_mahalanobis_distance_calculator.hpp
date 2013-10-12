#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_MAHALANOBIS_DISTANCE_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_MAHALANOBIS_DISTANCE_CALCULATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_y_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BAGCPConstrainedMahalanobisDistanceCalculator
{
public:
  typedef _Scalar Scalar;
  typedef BAGCPConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::YVector Residuals;
  typedef BAGCPConstrainedYCovarianceInverse<Scalar> YCovarianceInverse;

private:
  typedef typename VectorFunction::Index Index;
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef typename YCovarianceInverse::NaiveYCovarianceInverse::Block FeatureBlock;
  typedef typename YCovarianceInverse::GCPBlock GCPBlock;

public:
  Scalar operator() (const Residuals& residuals,
                     const YCovarianceInverse& y_covariance_inverse) const
  {
    Index number_of_features =
      y_covariance_inverse.naive_y_covariance_inverse.blocks.size();
    Scalar mahalanobis_distance = Scalar(0);
    for (Index i = 0; i < number_of_features; i++)
    {
      Vector2 feature = residuals.segment(
        i * VectorFunction::params_per_feature_,
        VectorFunction::params_per_feature_);
      const FeatureBlock& feature_block=
        y_covariance_inverse.naive_y_covariance_inverse.blocks[i];
      mahalanobis_distance += feature.transpose() * feature_block * feature;
    }
    Index feature_size =
      number_of_features * VectorFunction::params_per_feature_;
    Index number_of_gcps = Index(y_covariance_inverse.gcp_blocks.size());
    for (Index i = 0; i < number_of_gcps; i++)
    {
      Vector3 gcp = residuals.segment(
        feature_size + i * VectorFunction::params_per_point_,
        VectorFunction::params_per_point_);
      const GCPBlock& gcp_block = y_covariance_inverse.gcp_blocks[i];
      mahalanobis_distance += gcp.transpose() * gcp_block * gcp;
    }

    return mahalanobis_distance;
  }
};

}
}
}

#endif