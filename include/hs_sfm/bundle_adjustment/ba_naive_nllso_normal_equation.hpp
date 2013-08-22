#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NLLSO_NORMAL_EQUATION_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NLLSO_NORMAL_EQUATION_HPP_

#include <limits>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_math/linear_algebra/lafunc/fwd_decl.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar,
          typename _Index,
          int params_per_cam,
          int params_per_pt
          >
struct BANaiveNormalMatrix
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef _Index Index;
  typedef EIGEN_MAT(Scalar,
                    params_per_cam,
                    params_per_cam) CamBlock;
  typedef EIGEN_MAT(Scalar,
                    params_per_pt,
                    params_per_pt) PtBlock;
  typedef EIGEN_MAT(Scalar,
                    params_per_cam,
                    params_per_pt) MixBlock;

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
    Index number_of_camera_params = params_per_cam * number_of_cameras;
    Index number_of_point_params = params_per_pt * number_of_points;
    if (i >= number_of_camera_params + number_of_point_params ||
        j >= number_of_camera_params + number_of_point_params)
    {
      //不合法，返回0
      return Scalar(0);
    }
    else if (i < number_of_camera_params &&
             j < number_of_camera_params &&
             i / params_per_cam == j / params_per_cam)
    {
      //属于相机块
      Index camera_id = i / params_per_cam;
      Index camera_row_offset = i % params_per_cam;
      Index camera_col_offset = j % params_per_cam;
      return camera_blocks[camera_id](camera_row_offset, camera_col_offset);
    }
    else if (i >= number_of_camera_params &&
             j >= number_of_camera_params &&
             i / params_per_pt == j / params_per_pt)
    {
      //属于点块
      Index point_id = (i - number_of_camera_params) / params_per_pt;
      Index point_row_offset =
        (i - number_of_camera_params) % params_per_pt;
      Index point_col_offset =
        (j - number_of_camera_params) % params_per_pt;
      return point_blocks[point_id](point_row_offset, point_col_offset);
    }
    else if (i < number_of_camera_params &&
             j >= number_of_camera_params)
    {
      //属于混合块
      Index mix_row_id =
        i / params_per_cam;
      Index mix_col_id =
        (j - number_of_camera_params) / params_per_pt;
      Index mix_row_offset =
        i % params_per_cam;
      Index mix_col_offset =
        (j - number_of_camera_params) % params_per_pt;
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
        j / params_per_cam;
      Index mix_col_id =
        (i - number_of_camera_params) / params_per_pt;
      Index mix_row_offset =
        j % params_per_cam;
      Index mix_col_offset =
        (i - number_of_camera_params) % params_per_pt;
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

template <typename _Scalar,
          typename _Index,
          int params_per_cam,
          int params_per_pt>
struct BANaiveAugmentNormalMatrix
{
  typedef _Scalar Scalar;
  typedef _Index Index;
  typedef BANaiveNormalMatrix<Scalar, Index,
                              params_per_cam,
                              params_per_pt> NormalMatrix;

  BANaiveAugmentNormalMatrix(const NormalMatrix& N, Scalar mu)
    : ref_N(N), mu_(mu) {}

  const NormalMatrix& ref_N;
  Scalar mu_;
};

template <typename _Scalar,
          typename _Index,
          int params_per_cam,
          int params_per_pt>
struct BANaiveBVec
{
  typedef _Scalar Scalar;
  typedef _Index Index;
  typedef EIGEN_VEC(Scalar, params_per_cam) CamSegment;
  typedef EIGEN_VEC(Scalar, params_per_pt) PtSegment;
  typedef EIGEN_VECTOR(CamSegment) CamSegmentContainer;
  typedef EIGEN_VECTOR(PtSegment) PtSegmentContainer;

  Scalar operator[](Index i) const
  {
    Index camera_params_size = number_of_cameras * params_per_cam;
    Index params_size = camera_params_size + number_of_points * params_per_pt;
    if (i < 0 || i >= params_size)
    {
      return Scalar(0);
    }
    else if (i < camera_params_size)
    {
      return cam_segments[i / params_per_cam][i % params_per_cam];
    }
    else
    {
      Index off = i - camera_params_size;
      return pt_segments[off / params_per_pt][off % params_per_pt];
    }
  }

  Index number_of_cameras;
  Index number_of_points;
  CamSegmentContainer cam_segments;
  PtSegmentContainer pt_segments;
};

}
}

namespace math
{
namespace la
{

template <typename _Scalar,
          typename _Index,
          int params_per_cam,
          int params_per_pt>
class MatMaxDiagValFunc<hs::sfm::ba::BANaiveNormalMatrix<_Scalar, _Index,
                                                         params_per_cam,
                                                         params_per_pt> >
{
public:
  typedef _Scalar Scalar;
  typedef _Index Index;
  typedef hs::sfm::ba::BANaiveNormalMatrix<Scalar, Index,
                                           params_per_cam,
                                           params_per_pt> Mat;

  Scalar operator()(const Mat& N) const
  {
    Index params_size = params_per_cam * N.number_of_cameras +
                        params_per_pt * N.number_of_points;
    Scalar max_diag = std::numeric_limits<Scalar>::min();
    for (Index i = 0; i < params_size; i++)
    {
      max_diag = std::max(max_diag, N.coeff(i, i));
    }

    return max_diag;
  }
};

}
}
}

#endif
