#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_EQUATION_BUILDER_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_EQUATION_BUILDER_HPP_

#include <cassert>
#include <cmath>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_jac_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_feat_cov_inv.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_residuals.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
struct BANaiveNormalMatrix
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef BANaiveVecFunc<Scalar> VecFunc;
  typedef typename VecFunc::Index Index;
  typedef EIGEN_MAT(Scalar,
                    VecFunc::m_paramsPerCam,
                    VecFunc::m_paramsPerCam) CamBlock;
  typedef EIGEN_MAT(Scalar,
                    VecFunc::m_paramsPerPt,
                    VecFunc::m_paramsPerPt) PtBlock;
  typedef EIGEN_MAT(Scalar,
                    VecFunc::m_paramsPerCam,
                    VecFunc::m_paramsPerPt) MixBlock;

  typedef EIGEN_VECTOR(CamBlock) CamBlockContainer;
  typedef EIGEN_VECTOR(PtBlock) PtBlockContainer;
  typedef EIGEN_VECTOR(MixBlock) MixBlockContainer;

  enum
  {
    EigenDefaultMajor =
#if defined(EIGEN_DEFAULT_TO_ROW_MAJOR)
      Eigen::RowMajor
#else
      Eigen::ColMajor
#endif
  };
  typedef typename MixBlockContainer::size_type MixBlockIdx;
  typedef EIGEN_SMAT(MixBlockIdx, EigenDefaultMajor, Index) MixBlockMap;

  Scalar coeff(Index i, Index j) const
  {
    Index number_of_camera_params = VecFunc::m_paramsPerCam * number_of_cameras;
    Index number_of_point_params = VecFunc::m_paramsPerPt * number_of_points;
    if (i >= number_of_camera_params + number_of_point_params ||
        j >= number_of_camera_params + number_of_point_params)
    {
      //不合法，返回0
      return Scalar(0);
    }
    else if (i < number_of_camera_params &&
             j < number_of_camera_params &&
             i / VecFunc::m_paramsPerCam == j / VecFunc::m_paramsPerCam)
    {
      //属于相机块
      Index camera_id = i / VecFunc::m_paramsPerCam;
      Index camera_row_offset = i % VecFunc::m_paramsPerCam;
      Index camera_col_offset = j % VecFunc::m_paramsPerCam;
      return camera_blocks[camera_id](camera_row_offset, camera_col_offset);
    }
    else if (i >= number_of_camera_params &&
             j >= number_of_camera_params &&
             i / VecFunc::m_paramsPerPt == j / VecFunc::m_paramsPerPt)
    {
      //属于点块
      Index point_id = (i - number_of_camera_params) / VecFunc::m_paramsPerPt;
      Index point_row_offset =
        (i - number_of_camera_params) % VecFunc::m_paramsPerPt;
      Index point_col_offset =
        (j - number_of_camera_params) % VecFunc::m_paramsPerPt;
      return point_blocks[point_id](point_row_offset, point_col_offset);
    }
    else if (i < number_of_camera_params &&
             j >= number_of_camera_params)
    {
      //属于混合块
      Index mix_row_id =
        i / VecFunc::m_paramsPerCam;
      Index mix_col_id =
        (j - number_of_camera_params) / VecFunc::m_paramsPerPt;
      Index mix_row_offset =
        i % VecFunc::m_paramsPerCam;
      Index mix_col_offset =
        (j - number_of_camera_params) % VecFunc::m_paramsPerPt;
      MixBlockIdx mix_id =
        mix_block_map.coeff(mix_row_id, mix_col_id);
      if (mix_id == 0)
      {
        return Scalar(0);
      }
      else
      {
        return mix_blocks[mix_id - 1](mix_row_offset, mix_col_offset);
      }
    }
    else if (i >= number_of_camera_params &&
             j < number_of_camera_params)
    {
      //属于混合转置块
      Index mix_row_id =
        j / VecFunc::m_paramsPerCam;
      Index mix_col_id =
        (i - number_of_camera_params) / VecFunc::m_paramsPerPt;
      Index mix_row_offset =
        j % VecFunc::m_paramsPerCam;
      Index mix_col_offset =
        (i - number_of_camera_params) % VecFunc::m_paramsPerPt;
      MixBlockIdx mix_id =
        mix_block_map.coeff(mix_row_id, mix_col_id);
      if (mix_id == 0)
      {
        return Scalar(0);
      }
      else
      {
        return mix_blocks[mix_id - 1](mix_row_offset, mix_col_offset);
      }
    }
    else
    {
      //不属于任何块，为0
      return Scalar(0);
    }

  }

  Index number_of_cameras;
  Index number_of_points;
  CamBlockContainer camera_blocks;
  PtBlockContainer point_blocks;
  MixBlockContainer mix_blocks;
  MixBlockMap mix_block_map;
};

template <typename _Scalar>
struct BANaiveBVec
{
  typedef _Scalar Scalar;
  typedef BANaiveVecFunc<Scalar> VecFunc;
  typedef typename VecFunc::Index Index;
  typedef EIGEN_VEC(Scalar, VecFunc::m_paramsPerCam) CamSegment;
  typedef EIGEN_VEC(Scalar, VecFunc::m_paramsPerPt) PtSegment;
  typedef EIGEN_VECTOR(CamSegment) CamSegmentContainer;
  typedef EIGEN_VECTOR(PtSegment) PtSegmentContainer;

  Index number_of_cameras;
  Index number_of_points;
  CamSegmentContainer cam_segments;
  PtSegmentContainer pt_segments;
};

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

  typedef BANaiveNormalMatrix<Scalar> NormalMat;
  typedef BANaiveBVec<Scalar> BVec;

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
    //Camera blocks
    Index cam_num = j.m_camNum;
    Index pt_num = j.m_ptNum;
    N.camera_blocks.resize(cam_num, CamBlock::Zero());
    N.point_blocks.resize(pt_num, PtBlock::Zero());
    b.cam_segments.resize(cam_num, CamSegment::Zero());
    b.pt_segments.resize(pt_num, PtSegment::Zero());
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
      N.point_blocks[pt_drv_blk.m_ptId] +=
        pt_drv_blk.m_drvBlk.transpose() * y_cov_inv_blk *
        pt_drv_blk.m_drvBlk;
      b.pt_segments[pt_drv_blk.m_ptId] +=
        pt_drv_blk.m_drvBlk.transpose() * y_cov_inv_blk *
        r.segment(VecFunc::m_paramsPerFeat * featId,
                  VecFunc::m_paramsPerFeat);
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
};

}
}
}

#endif
