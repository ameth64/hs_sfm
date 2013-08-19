#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_RESIDUALS_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_RESIDUALS_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_feat_cov_inv.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VecFunc>
class BANaiveResidualsCalc;

template <typename _Scalar>
class BANaiveResidualsCalc<BANaiveVecFunc<_Scalar> >
{
public:
  typedef int Err;
  typedef _Scalar Scalar;
  typedef BANaiveVecFunc<Scalar> VecFunc;
  typedef typename VecFunc::XVec XVec;
  typedef typename VecFunc::YVec YVec;
  typedef BANaiveFeatCovInv<Scalar> YCovInv;
  typedef YVec Residuals;
  typedef typename VecFunc::Index Index;
  typedef EIGEN_VEC(Scalar, 2) Vec2;
  typedef typename YCovInv::Block YCovInvBlock;

  Err operator() (const VecFunc& f,
                  const XVec& x,
                  const YVec& y,
                  Residuals& r)
  {
    YVec noise_y;
    if (f(x, y) != 0)
    {
      return -1;
    }

    if (y.rows() != noise_y.rows())
    {
      return -1;
    }

    r = noise_y - y;

    return 0;
  }

  Residuals operator()(const YVec& noise_y, const YVec& estimated_y) const
  {
    return estimated_y - noise_y;
  }

  Scalar operator()(const Residuals& r, const YCovInv& y_cov_inv) const
  {
    Index feat_num = r.row() / VecFunc::m_paramsPerFeat;
    if (Index(y_cov_inv.m_blocks.size()) != feat_num)
    {
      //残差与协方差矩阵不匹配，返回0。
      return Scalar(0);
    }

    Scalar mahalanobis_distance = Scalar(0);
    for (Index i = 0; i < feat_num; i++)
    {
      Vec2 feat = r.segment(i * VecFunc::m_paramsPerFeat,
                            VecFunc::m_paramsPerFeat);

      YCovInvBlock y_cov_inv_block = y_cov_inv.m_blocks[i];
      mahalanobis_distance += feat.transpose() * y_cov_inv_block * feat;
    }

    return mahalanobis_distance;
  }
};

}
}
}

#endif
