#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_AUGMENTED_NORMAL_MATRIX_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_AUGMENTED_NORMAL_MATRIX_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_normal_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar,
          typename _Index,
          int params_per_camera,
          int params_per_point>
struct BANaiveAugmentedNormalMatrix
{
  typedef _Scalar Scalar;
  typedef _Index Index;
  typedef BANaiveNormalMatrix<Scalar, Index,
                              params_per_camera,
                              params_per_point> NormalMatrix;

  BANaiveAugmentedNormalMatrix(const NormalMatrix& normal_matrix, Scalar mu)
    : normal_matrix_ref_(normal_matrix), mu_(mu) {}

  const NormalMatrix& normal_matrix_ref_;
  Scalar mu_;
};

}
}
}

#endif