#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_FFD_JAC_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_FFD_JAC_HPP_

#include <cmath>

#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_jac_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VecFunc>
class BANaiveFFDJacobian;

/**
 *  利用Forward Finit Differace计算简单的Bundle Adjustment函数的Jacobian矩阵
 */
template <typename _Scalar>
class BANaiveFFDJacobian<BANaiveVecFunc<_Scalar> >
{
public:
  typedef _Scalar Scalar;

  typedef int Err;

  typedef BANaiveVecFunc<Scalar> VecFunc;
  typedef typename VecFunc::XVec XVec;
  typedef typename VecFunc::YVec YVec;
  typedef typename VecFunc::Index Index;
  typedef typename VecFunc::FeatMap FeatMap;
  typedef typename VecFunc::FeatMapContainer FeatMapContainer;
  typedef BANaiveJacMatrix<Scalar, Index,
                           VecFunc::m_paramsPerFeat,
                           VecFunc::m_paramsPerCam,
                           VecFunc::m_paramsPerPt> Jac;
  typedef typename Jac::DrvIdx DrvIdx;
  typedef EIGEN_VEC(Scalar, 3) Vec3;
  typedef EIGEN_VEC(Scalar, 2) Vec2;

  BANaiveFFDJacobian()
    : m_delta(Scalar(1e-4)), m_minDelta(Scalar(1e-6)) {}

  BANaiveFFDJacobian(Scalar delta, Scalar minDelta)
    : m_delta(delta), m_minDelta(minDelta) {}

  Err operator()(const VecFunc& vecFunc, const XVec& x, Jac& jac) const
  {
    typedef Eigen::Triplet<DrvIdx, Index> TripletType;
    std::vector<TripletType> camsTriplet;
    std::vector<TripletType> ptsTriplet;
    Index camNum = vecFunc.getCamNum();
    Index ptNum = vecFunc.getPtNum();
    Index featNum = vecFunc.getFeatNum();
    Index camParamsSize = vecFunc.getCamParamsSize();
    FeatMapContainer featMaps = vecFunc.getFeatMaps();
    Index xSize = vecFunc.getXSize();
    jac.clear();
    jac.m_camNum = camNum;
    jac.m_ptNum = ptNum;
    for (Index i = 0; i < featNum; i++)
    {
      const FeatMap& featMap = featMaps[i];
      Index j = featMap.first;
      Index k = featMap.second;

      Vec3 r = x.segment(j * VecFunc::m_paramsPerCam, 3);
      Vec3 t = x.segment(j * VecFunc::m_paramsPerCam + 3, 3);
      Vec3 p = x.segment(camParamsSize + k * VecFunc::m_paramsPerPt,
                         VecFunc::m_paramsPerPt);
      Vec2 feat;
      VecFunc::ptPrjToFeat(r, t, p, feat);

      Jac::CamDrvBlk camDrvBlk;
      camDrvBlk.m_camId = j;
      camDrvBlk.m_featId = i;
      Vec3 rDiff = r;
      for (Index m = 0; m < 3; m++)
      {
        Scalar d = std::max(std::abs(m_delta * r[m]),
                            m_minDelta);
        rDiff[m] += d;
        Vec2 featDiff;
        VecFunc::ptPrjToFeat(rDiff, t, p, featDiff);
        for (Index n = 0; n < VecFunc::m_paramsPerFeat; n++)
        {
          camDrvBlk.m_drvBlk.col(m)[n] = (featDiff[n] - feat[n]) / d;
        }
        rDiff[m] = r[m];
      }
      Vec3 tDiff = t;
      for (Index m = 0; m < 3; m++)
      {
        Scalar d = std::max(std::abs(m_delta * t[m]),
                            m_minDelta);
        tDiff[m] += d;
        Vec2 featDiff;
        VecFunc::ptPrjToFeat(r, tDiff, p, featDiff);
        for (Index n = 0; n < VecFunc::m_paramsPerFeat; n++)
        {
          camDrvBlk.m_drvBlk.col(3 + m)[n] = (featDiff[n] - feat[n]) / d;
        }
        tDiff[m] = t[m];
      }
      jac.m_camsDrv.push_back(camDrvBlk);
      camsTriplet.push_back(TripletType(j, k, jac.m_camsDrv.size()));

      Jac::PtDrvBlk ptDrvBlk;
      ptDrvBlk.m_ptId = k;
      ptDrvBlk.m_featId = i;
      Vec3 pDiff = p;
      for (Index m = 0; m < VecFunc::m_paramsPerPt; m++)
      {
        Scalar d = std::max(std::abs(m_delta * p[m]),
                            m_minDelta);
        pDiff[m] += d;
        Vec2 featDiff;
        VecFunc::ptPrjToFeat(r, t, pDiff, featDiff);
        for (Index n = 0; n < VecFunc::m_paramsPerFeat; n++)
        {
          ptDrvBlk.m_drvBlk.col(m)[n] = (featDiff[n] - feat[n]) / d;
        }
        pDiff[m] = p[m];
      }
      jac.m_ptsDrv.push_back(ptDrvBlk);
      ptsTriplet.push_back(TripletType(j, k, jac.m_ptsDrv.size()));
    }
    jac.m_camsDrvMap.resize(camNum, ptNum);
    jac.m_camsDrvMap.setFromTriplets(camsTriplet.begin(), 
                     camsTriplet.end());
    jac.m_ptsDrvMap.resize(camNum, ptNum);
    jac.m_ptsDrvMap.setFromTriplets(ptsTriplet.begin(),
                    ptsTriplet.end());

    return 0;
  }

private:
  Scalar m_delta;
  Scalar m_minDelta;
};

}
}
}

#endif
