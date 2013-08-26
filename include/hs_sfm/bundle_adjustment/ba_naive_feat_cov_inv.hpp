#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_FEAT_COV_INV_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_FEAT_COV_INV_HPP_

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
class BANaiveFeatCovInv
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef EIGEN_MAT(Scalar, 2, 2) Block;
  typedef EIGEN_VECTOR(Block) BlockContainer;

  BlockContainer m_blocks;
};

}
}
}

#endif
