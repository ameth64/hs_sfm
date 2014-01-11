#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_MATRIX_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_MATRIX_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar,
          typename _Index,
          int params_per_camera,
          int params_per_point
          >
struct BANaiveNormalMatrix
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef _Index Index;
  typedef EIGEN_MATRIX(Scalar,
                       params_per_camera,
                       params_per_camera) CameraBlock;
  typedef EIGEN_MATRIX(Scalar,
                       params_per_point,
                       params_per_point) PointBlock;
  typedef EIGEN_MATRIX(Scalar,
                       params_per_camera,
                       params_per_point) MixBlock;

  typedef EIGEN_STD_VECTOR(CameraBlock) CameraBlockContainer;
  typedef EIGEN_STD_VECTOR(PointBlock) PointBlockContainer;
  typedef EIGEN_STD_VECTOR(MixBlock) MixBlockContainer;

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
  typedef EIGEN_SPARSE_MATRIX(MixBlockIdx, EigenDefaultMajor, Index)
          MixBlockMap;

  Scalar coeff(Index i, Index j) const
  {
    Index number_of_camera_params = params_per_camera * number_of_cameras;
    Index number_of_point_params = params_per_point * number_of_points;
    if (i >= number_of_camera_params + number_of_point_params ||
        j >= number_of_camera_params + number_of_point_params)
    {
      //不合法，返回0
      return Scalar(0);
    }
    else if (i < number_of_camera_params &&
             j < number_of_camera_params &&
             i / params_per_camera == j / params_per_camera)
    {
      //属于相机块
      Index camera_id = i / params_per_camera;
      Index camera_row_offset = i % params_per_camera;
      Index camera_col_offset = j % params_per_camera;
      return camera_blocks[camera_id](camera_row_offset, camera_col_offset);
    }
    else if (i >= number_of_camera_params &&
             j >= number_of_camera_params &&
             i / params_per_point == j / params_per_point)
    {
      //属于点块
      Index point_id = (i - number_of_camera_params) / params_per_point;
      Index point_row_offset =
        (i - number_of_camera_params) % params_per_point;
      Index point_col_offset =
        (j - number_of_camera_params) % params_per_point;
      return point_blocks[point_id](point_row_offset, point_col_offset);
    }
    else if (i < number_of_camera_params &&
             j >= number_of_camera_params)
    {
      //属于混合块
      Index mix_row_id =
        i / params_per_camera;
      Index mix_col_id =
        (j - number_of_camera_params) / params_per_point;
      Index mix_row_offset =
        i % params_per_camera;
      Index mix_col_offset =
        (j - number_of_camera_params) % params_per_point;
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
        j / params_per_camera;
      Index mix_col_id =
        (i - number_of_camera_params) / params_per_point;
      Index mix_row_offset =
        j % params_per_camera;
      Index mix_col_offset =
        (i - number_of_camera_params) % params_per_point;
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
  CameraBlockContainer camera_blocks;
  PointBlockContainer point_blocks;
  MixBlockContainer mix_blocks;
  MixBlockMap mix_block_map;
};

}
}
}

#endif
