#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_AUGMENTED_NORMAL_MATRIX_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_AUGMENTED_NORMAL_MATRIX_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_normal_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

/**
 *  多个相机共享的Bundle Adjustment向量函数的增广正规矩阵。
 */
template <typename _Scalar>
class CameraSharedAugmentedNormalMatrix
{
public:
  typedef _Scalar Scalar;
  typedef CameraSharedNormalMatrix<Scalar> NormalMatrix;

  CameraSharedAugmentedNormalMatrix(const NormalMatrix& normal_matrix,
                                    Scalar mu)
    : normal_matrix_ref_(normal_matrix), mu_(mu) {}

  const NormalMatrix& normal_matrix_ref() const
  {
    return normal_matrix_ref_;
  }

  Scalar mu() const
  {
    return mu_;
  }

private:
  const NormalMatrix& normal_matrix_ref_;
  Scalar mu_;
};

}
}
}

#endif
