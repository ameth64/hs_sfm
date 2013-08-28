#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_EQUATION_BUILDER_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_EQUATION_BUILDER_HPP_

#include <cassert>
#include <cmath>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_feat_cov_inv.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_jac_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_nllso_normal_equation.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_residuals.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BANaiveNormalEquationBuilder
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef BANaiveVecFunc<Scalar> VecFunc;
  typedef typename VecFunc::XVec XVec;
  typedef typename VecFunc::YVec YVec;
  typedef typename VecFunc::Index Index;
  typedef BANaiveJacMatrix<Scalar, Index,
          VecFunc::m_paramsPerFeat,
          VecFunc::m_paramsPerCam,
          VecFunc::m_paramsPerPt> Jac;
  typedef BANaiveFeatCovInv<Scalar> YCovInv;

  typedef BANaiveResidualsCalc<VecFunc> ResidualsCalc;
  typedef typename ResidualsCalc::Residuals Residuals;

  typedef BANaiveNormalMatrix<Scalar, Index,
                              VecFunc::m_paramsPerCam,
                              VecFunc::m_paramsPerPt> NormalMat;
  typedef BANaiveBVec<Scalar, Index,
                      VecFunc::m_paramsPerCam,
                      VecFunc::m_paramsPerPt> BVec;

  typedef typename Jac::CamDrvBlk CamDrvBlk;
  typedef typename Jac::PtDrvBlk PtDrvBlk;
  typedef typename Jac::DrvIdx DrvIdx;
  typedef typename Jac::DerivativeMap DrivativeMap;
  typedef typename YCovInv::Block YCovInvBlk;
  typedef typename NormalMat::CamBlock CamBlock;
  typedef typename NormalMat::PtBlock PtBlock;
  typedef typename NormalMat::MixBlock MixBlock;
  typedef typename NormalMat::MixBlockIdx MixBlockIdx;
  typedef typename NormalMat::MixBlockMap MixBlockMap;
  typedef typename BVec::CamSegment CamSegment;
  typedef typename BVec::PtSegment PtSegment;

  Err operator()(const Jac& j, const Residuals& r, const YCovInv& ycovInv,
                 NormalMat& N, BVec& b) const
  {
    if (BuildNormalMatrix(j, ycovInv, N) != 0) return -1;
    if (BuildBVector(j, r, ycovInv, b) != 0) return -1;

    return 0;
  }

  Err BuildNormalMatrix(const Jac& j, const YCovInv& ycovInv,
                        NormalMat& N) const
  {
    //Camera blocks
    Index cam_num = j.m_camNum;
    Index pt_num = j.m_ptNum;
    N.camera_blocks.resize(cam_num, CamBlock::Zero());
    N.point_blocks.resize(pt_num, PtBlock::Zero());
    size_t cam_blks_num = j.m_camsDrv.size();
    for (size_t i = 0; i < cam_blks_num; i++)
    {
      const CamDrvBlk& cam_drv_blk = j.m_camsDrv[i];
      Index featId = cam_drv_blk.m_featId;
      const YCovInvBlk& y_cov_inv_blk =
        ycovInv.m_blocks[featId];
      N.camera_blocks[cam_drv_blk.m_camId] +=
        cam_drv_blk.m_drvBlk.transpose() * y_cov_inv_blk *
        cam_drv_blk.m_drvBlk;
    }

    size_t pt_blks_num = j.m_ptsDrv.size();
    for (size_t i = 0; i < pt_blks_num; i++)
    {
      const PtDrvBlk& pt_drv_blk = j.m_ptsDrv[i];
      Index featId = pt_drv_blk.m_featId;
      const YCovInvBlk& y_cov_inv_blk =
        ycovInv.m_blocks[featId];
      N.point_blocks[pt_drv_blk.m_ptId] +=
        pt_drv_blk.m_drvBlk.transpose() * y_cov_inv_blk *
        pt_drv_blk.m_drvBlk;
    }

    typedef Eigen::Triplet<MixBlockIdx, Index> TripletType;
    std::vector<TripletType> mix_block_triplets;
    N.mix_blocks.clear();
    N.number_of_cameras = cam_num;
    N.number_of_points = pt_num;
    for (Index i = 0; i < cam_num; i++)
    {
      for (Index k = 0; k < pt_num; k++)
      {
        DrvIdx cam_id = j.m_camsDrvMap.coeff(i, k);
        DrvIdx pt_id = j.m_ptsDrvMap.coeff(i, k);
        if (cam_id != 0 &&
            pt_id != 0)
        {
          const CamDrvBlk& cam_drv_blk = j.m_camsDrv[cam_id - 1];
          const PtDrvBlk& pt_drv_blk = j.m_ptsDrv[pt_id - 1];
          assert(cam_drv_blk.m_featId == pt_drv_blk.m_featId);
          Index featId = cam_drv_blk.m_featId;
          const YCovInvBlk& y_cov_inv_blk =
            ycovInv.m_blocks[featId];
          MixBlock mix_block =
            cam_drv_blk.m_drvBlk.transpose() *
            y_cov_inv_blk *
            pt_drv_blk.m_drvBlk;
          N.mix_blocks.push_back(mix_block);
          mix_block_triplets.push_back(TripletType(i, k, N.mix_blocks.size()));
        }
      }
    }
    N.mix_block_map.resize(cam_num, pt_num);
    N.mix_block_map.setFromTriplets(mix_block_triplets.begin(),
                                    mix_block_triplets.end());

    return 0;
  }

  Err BuildBVector(const Jac& j, const Residuals& r, const YCovInv& ycovInv,
                   BVec& b) const
  {
    //Camera blocks
    Index cam_num = j.m_camNum;
    Index pt_num = j.m_ptNum;
    b.cam_segments.resize(cam_num, CamSegment::Zero());
    b.pt_segments.resize(pt_num, PtSegment::Zero());
    size_t cam_blks_num = j.m_camsDrv.size();
    for (size_t i = 0; i < cam_blks_num; i++)
    {
      const CamDrvBlk& cam_drv_blk = j.m_camsDrv[i];
      Index featId = cam_drv_blk.m_featId;
      const YCovInvBlk& y_cov_inv_blk =
        ycovInv.m_blocks[featId];
      b.cam_segments[cam_drv_blk.m_camId] +=
        cam_drv_blk.m_drvBlk.transpose() * y_cov_inv_blk *
        r.segment(VecFunc::m_paramsPerFeat * featId,
                  VecFunc::m_paramsPerFeat);
    }

    size_t pt_blks_num = j.m_ptsDrv.size();
    for (size_t i = 0; i < pt_blks_num; i++)
    {
      const PtDrvBlk& pt_drv_blk = j.m_ptsDrv[i];
      Index featId = pt_drv_blk.m_featId;
      const YCovInvBlk& y_cov_inv_blk =
        ycovInv.m_blocks[featId];
      b.pt_segments[pt_drv_blk.m_ptId] +=
        pt_drv_blk.m_drvBlk.transpose() * y_cov_inv_blk *
        r.segment(VecFunc::m_paramsPerFeat * featId,
                  VecFunc::m_paramsPerFeat);
    }

    b.number_of_cameras = cam_num;
    b.number_of_points = pt_num;

    return 0;
  }
};

}
}
}

#endif
