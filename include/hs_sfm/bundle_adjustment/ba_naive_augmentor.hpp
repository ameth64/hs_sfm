#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_AUGMENTOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_AUGMENTOR_HPP_

#include "hs_sfm/bundle_adjustment/ba_naive_augmented_normal_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BANaiveAugmentor
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef BANaiveAugmentedNormalMatrix<Scalar, Index,
                                       VectorFunction::params_per_camera_,
                                       VectorFunction::params_per_point_>
          AugmentedNormalMatrix;
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
