#ifndef _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_Y_COVARIANCE_INVERSE_HPP_
#define _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_Y_COVARIANCE_INVERSE_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

template <typename _Scalar>
struct IntrinsicConstrainedYCovarianceInverse
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef EIGEN_MATRIX(Scalar, 2, 2) KeyBlock;
  typedef EIGEN_STD_VECTOR(KeyBlock) KeyBlockContainer;

  KeyBlockContainer key_blocks;
  Scalar focal_length_stddev;
  Scalar skew_stddev;
  Scalar principal_point_x_stddev;
  Scalar principal_point_y_stddev;
  Scalar pixel_ratio_stddev;
  Scalar k1_stddev;
  Scalar k2_stddev;
  Scalar k3_stddev;
  Scalar d1_stddev;
  Scalar d2_stddev;
};

}
}
}

#endif