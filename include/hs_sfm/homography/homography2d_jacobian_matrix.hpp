#ifndef _HS_SFM_HOMOGRAPHY_HOMOGRAPHY2D_JACOBIAN_MATRIX_HPP_
#define _HS_SFM_HOMOGRAPHY_HOMOGRAPHY2D_JACOBIAN_MATRIX_HPP_

#include <limits>
#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace homography
{

template <typename _Scalar>
class Homography2DJacobianMatrix
{
public:
  typedef _Scalar Scalar;
  typedef EIGEN_MATRIX(Scalar, 2, 2) KeyBlock;
  typedef EIGEN_STD_VECTOR(KeyBlock) KeyBlockContainer;
  typedef EIGEN_MATRIX(Scalar, 2, 9) HBlock;
  typedef EIGEN_STD_VECTOR(HBlock) HBlockContainer;
  typedef typename HBlock::Index Index;

public:
  const KeyBlockContainer& key_blocks() const {return key_blocks_;}
  KeyBlockContainer& key_blocks() {return key_blocks_;}
  const HBlockContainer& h_blocks() const {return h_blocks_;}
  HBlockContainer& h_blocks() {return h_blocks_;}

  Scalar coeff(Index i, Index j) const
  {
    Index key_size = Index(key_blocks_.size());
    Index h_size = Index(h_blocks_.size());
    if (key_size != h_size)
    {
      return std::numeric_limits<Scalar>::quiet_NaN();
    }
    if (j >= key_size * 2)
    {
      if (i % 4 > 1)
      {
        Index key_id = i / 4;
        Index key_offset = (i % 4) - 2;
        Index h_offset = j - key_size * 2;
        return h_blocks_[key_id](key_offset, h_offset);
      }
      else
      {
        return Scalar(0);
      }
    }
    else
    {
      Index x_key_id = j / 2;
      Index x_key_offset = j % 2;
      Index y_key_id = i / 4;
      if (x_key_id == y_key_id)
      {
        if (i % 4 < 2)
        {
          return ((i % 4) == x_key_offset ? Scalar(1) : Scalar(0));
        }
        else
        {
          Index y_key_offset = (i % 4) - 2;
          return key_blocks_[x_key_id](y_key_offset, x_key_offset);
        }
      }
      else
      {
        return Scalar(0);
      }
    }
  }

private:
  KeyBlockContainer key_blocks_;
  HBlockContainer h_blocks_;
};

}
}
}

#endif