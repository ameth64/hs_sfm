#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_FEATURE_COVARIANCE_INVERSE_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_FEATURE_COVARIANCE_INVERSE_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

/**
 *  简单BA函数的特征点协方差矩阵的逆矩阵。
 *  TODO:目前没提供接口，直接使用里边的数据，应该封装一下。
 */
template <typename _Scalar>
struct BANaiveFeatureCovarianceInverse
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef EIGEN_MATRIX(Scalar, 2, 2) Block;
  typedef EIGEN_STD_VECTOR(Block) BlockContainer;

  BlockContainer blocks;
};

}
}
}

#endif
