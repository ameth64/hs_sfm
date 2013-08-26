#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NLLSO_AUGMENTOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NLLSO_AUGMENTOR_HPP_

#include "hs_sfm/bundle_adjustment/ba_naive_nllso_normal_equation.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"

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
  typedef BANaiveVecFunc<Scalar> VecFunc;
  typedef typename VecFunc::Index Index;
  typedef BANaiveAugmentNormalMatrix<Scalar, Index,
                                     VecFunc::m_paramsPerCam,
                                     VecFunc::m_paramsPerPt> AugNormalMat;
  typedef typename AugNormalMat::NormalMatrix NormalMat;

  AugNormalMat operator()(const NormalMat& N, Scalar mu) const
  {
    return AugNormalMat(N, mu);
  }
};

}
}
}

#endif
