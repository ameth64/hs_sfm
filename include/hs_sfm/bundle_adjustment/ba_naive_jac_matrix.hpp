#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_JAC_MATRIX_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_JAC_MATRIX_HPP_

#include <vector>

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar,
          typename _Index,
          _Index paramsPerFeat,
          _Index paramsPerCam,
          _Index paramsPerPt>
struct BANaiveJacMatrix
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef _Index Index;


  struct CamDrvBlk
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef EIGEN_MAT(Scalar, paramsPerFeat,
                  paramsPerCam) DrvBlk;
    Index m_camId;
    Index m_featId;
    DrvBlk m_drvBlk;
  };

  struct PtDrvBlk
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef EIGEN_MAT(Scalar, paramsPerFeat,
                  paramsPerPt) DrvBlk;
    Index m_ptId;
    Index m_featId;
    DrvBlk m_drvBlk;
  };

  typedef EIGEN_VECTOR(CamDrvBlk) CamsDrv;
  typedef EIGEN_VECTOR(PtDrvBlk) PtsDrv;
  enum
  {
    EigenDefaultMajor =
#if defined(EIGEN_DEFAULT_TO_ROW_MAJOR)
    Eigen::RowMajor
#else
    Eigen::ColMajor
#endif
  };
  typedef typename CamsDrv::size_type DrvIdx;
  typedef EIGEN_SMAT(DrvIdx, EigenDefaultMajor, Index) DerivativeMap;

  inline void clear()
  {
    m_camsDrv.clear();
    m_camsDrvMap.swap(DerivativeMap());
    m_ptsDrv.clear();
    m_ptsDrvMap.swap(DerivativeMap());
  }

  Scalar coeff(Index i, Index j) const
  {
    Index off = j - paramsPerCam * m_camNum;
    Index featId = i / paramsPerFeat;
    Index featParamId = i % paramsPerFeat;
    if (off < 0)
    {
      //cam param
      Index camId = j / paramsPerCam;
      Index camParamId = j % paramsPerCam;
      const CamDrvBlk& camDrvBlk = m_camsDrv[featId];
      if (camDrvBlk.m_camId == camId)
      {
        return camDrvBlk.m_drvBlk(featParamId, camParamId);
      }
      else
      {
        return Scalar(0);
      }
    }
    else
    {
      //pt param
      Index ptId = off / paramsPerPt;
      Index ptParamId = off % paramsPerPt;
      const PtDrvBlk& ptDrvBlk = m_ptsDrv[featId];
      if (ptDrvBlk.m_ptId == ptId)
      {
        return ptDrvBlk.m_drvBlk(featParamId, ptParamId);
      }
      else
      {
        return Scalar(0);
      }
    }
  }

  CamsDrv m_camsDrv;
  DerivativeMap m_camsDrvMap;
  PtsDrv m_ptsDrv;
  DerivativeMap m_ptsDrvMap;
  Index m_camNum;
  Index m_ptNum;
};

}
}
}

#endif
