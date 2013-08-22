#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NLLSO_NORMAL_EQUATION_SOLVER_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NLLSO_NORMAL_EQUATION_SOLVER_HPP_

#include "hs_sfm/bundle_adjustment/ba_naive_nllso_normal_equation.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BANaiveNormalEquationSolver
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef BANaiveVecFunc<Scalar> VecFunc;
  typedef typename VecFunc::XVec XVec;
  typedef typename VecFunc::YVec YVec;
  typedef typename VecFunc::Index Index;
  typedef BANaiveNormalMatrix<Scalar, Index,
                              VecFunc::m_paramsPerCam,
                              VecFunc::m_paramsPerPt> NormalMat;
  typedef BANaiveBVec<Scalar, Index,
                      VecFunc::m_paramsPerCam,
                      VecFunc::m_paramsPerPt> BVec;

  typedef typename NormalMat::CamBlock CamBlock;
  typedef typename NormalMat::PtBlock PtBlock;
  typedef typename NormalMat::MixBlock MixBlock;
  typedef typename NormalMat::MixBlockIdx MixBlockIdx;
  typedef typename NormalMat::MixBlockMap MixBlockMap;
  typedef typename BVec::CamSegment CamSegment;
  typedef typename BVec::PtSegment PtSegment;

  enum
  {
    EigenDefaultMajor =
#if defined(EIGEN_DEFAULT_TO_ROW_MAJOR)
      Eigen::RowMajor
#else
      Eigen::ColMajor
#endif
  };
  typedef EIGEN_SMAT(Scalar, EigenDefaultMajor, Index)
          SchurMat;
  typedef EIGEN_VEC(Scalar, Eigen::Dynamic) CameraParams;
  typedef EIGEN_VEC(Scalar, Eigen::Dynamic) PointParams;
  typedef EIGEN_VEC(Scalar, Eigen::Dynamic) CameraLHS;
  typedef EIGEN_MAT(Scalar,
                    VecFunc::m_paramsPerCam,
                    VecFunc::m_paramsPerPt) MixPointBlock;
  typedef EIGEN_VECTOR(MixPointBlock) MixPointBlockContainer;
  typedef typename MixPointBlockContainer::size_type MixPointBlockIdx;
  typedef EIGEN_SMAT(MixPointBlockIdx, EigenDefaultMajor, Index) MixPointBlockMap;

  Err operator()(const NormalMat& N, const BVec& b, XVec& x) const
  {
    if (SolveCameraParams(N, b, x) != 0) return -1;
    if (SolvePointParams(N, b, x) != 0) return -1;
    return 0;
  }

private:
  Err SolveCameraParams(const NormalMat& N, const BVec& b, XVec& x) const
  {
    Index schur_mat_size = N.number_of_cameras * VecFunc::m_paramsPerCam;
    SchurMat schur_mat;
    x.resize(schur_mat_size + N.number_of_points * VecFunc::m_paramsPerPt);

    //compute mix point blocks
    MixPointBlockContainer mix_point_blocks;
    MixPointBlockMap mix_point_block_map;
    typedef Eigen::Triplet<MixPointBlockIdx, Index> MixPointTripletType;
    std::vector<MixPointTripletType> mix_point_block_triplets;
    for (Index i = 0; i < N.number_of_cameras; i++)
    {
      for (Index j = 0; j < N.number_of_points; j++)
      {
        MixBlockIdx mix_id = N.mix_block_map.coeff(i, j);
        if (mix_id > 0)
        {
          const MixBlock& mix = N.mix_blocks[mix_id - 1];
          const PtBlock& point = N.point_blocks[j];
          MixPointBlock mix_point = mix * point.inverse();
          mix_point_blocks.push_back(mix_point);
          mix_point_block_triplets.push_back(
            MixPointTripletType(i, j, mix_point_blocks.size()));
        }
      }
    }
    mix_point_block_map.resize(N.number_of_cameras,
                               N.number_of_points);
    mix_point_block_map.setFromTriplets(mix_point_block_triplets.begin(),
                                        mix_point_block_triplets.end());

    //compute schur matrix
    typedef Eigen::Triplet<Scalar, Index> SchurTripletType;
    std::vector<SchurTripletType> schur_triplets;
    for (Index i = 0; i < N.number_of_cameras; i++)
    {
      for (Index j = i; j < N.number_of_cameras; j++)
      {
        CamBlock camera_block = CamBlock::Zero();
        for (Index k = 0; k < N.number_of_points; k++)
        {
          MixPointBlockIdx mix_point_id = mix_point_block_map.coeff(i, k);
          MixBlockIdx mix_id = N.mix_block_map.coeff(j, k);
          if (mix_point_id != 0 &&
              mix_id != 0)
          {
            camera_block -= mix_point_blocks[mix_point_id - 1] *
                            N.mix_blocks[mix_id - 1].transpose();
          }
        }

        //fill schur matrix
        if (i == j)
        {
          camera_block += N.camera_blocks[i];
          for (Index m = 0; m < VecFunc::m_paramsPerCam; m++)
          {
            for (Index n = 0; n < VecFunc::m_paramsPerCam; n++)
            {
              schur_triplets.push_back(SchurTripletType(
                                        i * VecFunc::m_paramsPerCam + m,
                                        j * VecFunc::m_paramsPerCam + n,
                                        camera_block(m, n)));
            }
          }
        }
        else
        {
          for (Index m = 0; m < VecFunc::m_paramsPerCam; m++)
          {
            for (Index n = 0; n < VecFunc::m_paramsPerCam; n++)
            {
              schur_triplets.push_back(SchurTripletType(
                                        i * VecFunc::m_paramsPerCam + m,
                                        j * VecFunc::m_paramsPerCam + n,
                                        camera_block(m, n)));
              schur_triplets.push_back(SchurTripletType(
                                        j * VecFunc::m_paramsPerCam + n,
                                        i * VecFunc::m_paramsPerCam + m,
                                        camera_block(m, n)));
            }
          }
        }
      }
    }
    schur_mat.resize(schur_mat_size, schur_mat_size);
    schur_mat.setFromTriplets(schur_triplets.begin(),
                              schur_triplets.end());

    //compute camera LHS
    CameraLHS camera_lhs(schur_mat_size);
    for (Index i = 0; i < N.number_of_cameras; i++)
    {
      CamSegment cam_segment = b.cam_segments[i];
      for (Index j = 0; j < N.number_of_points; j++)
      {
        MixPointBlockIdx mix_point_id = mix_point_block_map.coeff(i, j);
        if (mix_point_id > 0)
        {
          const PtSegment& point_segment = b.pt_segments[j];
          cam_segment -= mix_point_blocks[mix_point_id - 1] *
                         point_segment;
        }
      }
      camera_lhs.segment(i * VecFunc::m_paramsPerCam,
                         VecFunc::m_paramsPerCam) = cam_segment;
    }

    //solve camera params
    typedef Eigen::SimplicialLDLT<SchurMat> Solver;
    Solver solver;
    solver.compute(schur_mat);
    if (solver.info() != Eigen::Success)
    {
      return -1;
    }
    x.segment(0, schur_mat_size) = solver.solve(camera_lhs);
    if (solver.info() != Eigen::Success)
    {
      return -1;
    }

    return 0;
  }

  Err SolvePointParams(const NormalMat& N, const BVec& b, XVec& x) const
  {
    Index camera_params_size = N.number_of_cameras * VecFunc::m_paramsPerCam;
    for (Index i = 0; i < N.number_of_points; i++)
    {
      PtSegment pt_segment = b.pt_segments[i];
      for (Index j = 0; j < N.number_of_cameras; j++)
      {
        MixBlockIdx mix_id = N.mix_block_map.coeff(j, i);
        if (mix_id > 0)
        {
          const MixBlock& mix = N.mix_blocks[mix_id - 1];
          pt_segment -= mix.transpose() *
                        x.segment(j * VecFunc::m_paramsPerCam,
                                  VecFunc::m_paramsPerCam);
        }
      }
      x.segment(camera_params_size + i * VecFunc::m_paramsPerPt,
                VecFunc::m_paramsPerPt) =
        N.point_blocks[i].inverse() * pt_segment;
    }

    return 0;
  }
};

}
}
}

#endif
