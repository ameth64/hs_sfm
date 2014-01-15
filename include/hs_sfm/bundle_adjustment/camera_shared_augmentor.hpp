#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_AUGMENTOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_AUGMENTOR_HPP_

#include "hs_sfm/bundle_adjustment/camera_shared_augmented_normal_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedAugmentor
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef CameraSharedAugmentedNormalMatrix<Scalar> AugmentedNormalMatrix;
  typedef typename AugmentedNormalMatrix::NormalMatrix NormalMatrix;

  AugmentedNormalMatrix operator()(const NormalMatrix& normal_matrix,
                                   Scalar mu) const
  {
    return AugmentedNormalMatrix(normal_matrix, mu);
  }
};

}
}
}

#endif
