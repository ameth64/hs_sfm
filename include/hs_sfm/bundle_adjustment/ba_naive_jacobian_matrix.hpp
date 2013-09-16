#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_JACOBIAN_MATRIX_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_JACOBIAN_MATRIX_HPP_

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
          _Index params_per_feature,
          _Index params_per_camera,
          _Index params_per_point>
class BANaiveJacobianMatrix
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef _Index Index;


  struct CameraDerivativeBlock
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef EIGEN_MATRIX(Scalar, params_per_feature,  params_per_camera)
            DerivativeBlock;
    Index camera_id;
    Index feature_id;
    DerivativeBlock derivative_block;
  };

  struct PointDerivativeBlock
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef EIGEN_MATRIX(Scalar, params_per_feature, params_per_point)
            DerivativeBlock;
    Index point_id;
    Index feature_id;
    DerivativeBlock derivative_block;
  };

  typedef EIGEN_STD_VECTOR(CameraDerivativeBlock) CameraDerivativeContainer;
  typedef EIGEN_STD_VECTOR(PointDerivativeBlock) PointDerivativeContainer;
  enum
  {
    EigenDefaultMajor =
#if defined(EIGEN_DEFAULT_TO_ROW_MAJOR)
    Eigen::RowMajor
#else
    Eigen::ColMajor
#endif
  };
  typedef typename CameraDerivativeContainer::size_type DerivativeId;
  typedef EIGEN_SPARSE_MATRIX(DerivativeId, EigenDefaultMajor, Index)
          DerivativeMap;

  inline void clear()
  {
    CameraDerivativeContainer().swap(camera_derivatives_);
    PointDerivativeContainer().swap(point_derivatives_);
    DerivativeMap().swap(camera_derivatives_map_);
    DerivativeMap().swap(points_derivatives_map_);
  }

  Scalar coeff(Index i, Index j) const
  {
    Index off = j - params_per_camera * number_of_cameras_;
    Index feature_id = i / params_per_feature;
    Index feature_param_id = i % params_per_feature;
    if (off < 0)
    {
      //cam param
      Index camera_id = j / params_per_camera;
      Index camera_param_id = j % params_per_camera;
      const CameraDerivativeBlock& camera_derivative_block =
        camera_derivatives_[feature_id];
      if (camera_derivative_block.camera_id == camera_id)
      {
        return camera_derivative_block.derivative_block(
                 feature_param_id, camera_param_id);
      }
      else
      {
        return Scalar(0);
      }
    }
    else
    {
      //pt param
      Index point_id = off / params_per_point;
      Index point_param_id = off % params_per_point;
      const PointDerivativeBlock& point_derivative_block =
        point_derivatives_[feature_id];
      if (point_derivative_block.point_id == point_id)
      {
        return point_derivative_block.derivative_block(
          feature_param_id, point_param_id);
      }
      else
      {
        return Scalar(0);
      }
    }
  }

  Index number_of_cameras() const
  {
    return number_of_cameras_;
  }

  void set_number_of_cameras(Index number_of_cameras)
  {
    number_of_cameras_ = number_of_cameras;
  }

  Index number_of_points() const
  {
    return number_of_points_;
  }

  void set_number_of_points(Index number_of_points)
  {
    number_of_points_ = number_of_points;
  }

  const CameraDerivativeContainer& camera_derivatives() const
  {
    return camera_derivatives_;
  }

  CameraDerivativeContainer& camera_derivatives()
  {
    return camera_derivatives_;
  }

  const PointDerivativeContainer& point_derivatives() const
  {
    return point_derivatives_;
  }

  PointDerivativeContainer& point_derivatives()
  {
    return point_derivatives_;
  }

  const DerivativeMap& camera_derivatives_map() const
  {
    return camera_derivatives_map_;
  }

  DerivativeMap& camera_derivatives_map()
  {
    return camera_derivatives_map_;
  }

  const DerivativeMap& point_derivatives_map() const
  {
    return points_derivatives_map_;
  }

  DerivativeMap& point_derivatives_map()
  {
    return points_derivatives_map_;
  }

private:
  CameraDerivativeContainer camera_derivatives_;
  DerivativeMap camera_derivatives_map_;
  PointDerivativeContainer point_derivatives_;
  DerivativeMap points_derivatives_map_;
  Index number_of_cameras_;
  Index number_of_points_;
};

}
}
}

#endif
