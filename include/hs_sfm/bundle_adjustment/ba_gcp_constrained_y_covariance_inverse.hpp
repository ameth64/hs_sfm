#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_Y_COVARIANCE_INVERSE_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_Y_COVARIANCE_INVERSE_HPP_

#include "hs_sfm/bundle_adjustment/ba_naive_feature_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
struct BAGCPConstrainedYCovarianceInverse
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef BANaiveFeatureCovarianceInverse<Scalar> NaiveYCovarianceInverse;
  typedef EIGEN_MATRIX(Scalar, 3, 3) GCPBlock;
  typedef EIGEN_STD_VECTOR(GCPBlock) GCPBlockContainers;

  GCPBlockContainers gcp_blocks;
  NaiveYCovarianceInverse naive_y_covariance_inverse;
};

}
}
}

#endif