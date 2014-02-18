#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NORMAL_MATRIX_MAX_DIAGONAL_VALUE_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NORMAL_MATRIX_MAX_DIAGONAL_VALUE_CALCULATOR_HPP_

#include <algorithm>
#include <limits>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_normal_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedNormalMatrixMaxDiagonalValueCalculator
{
public:
  typedef _Scalar Scalar;
  typedef CameraSharedNormalMatrix<Scalar> NormalMatrix;
  typedef typename NormalMatrix::Index Index;

  Scalar operator() (const NormalMatrix& normal_matrix) const
  {
    Index x_size = normal_matrix.GetXSize();
    Scalar max_value = -std::numeric_limits<Scalar>::max();
    for (Index i = 0; i < x_size; i++)
    {
      max_value = std::max(max_value, normal_matrix.coeff(i, i));
    }

    return max_value;
  }
};

}
}
}

#endif
