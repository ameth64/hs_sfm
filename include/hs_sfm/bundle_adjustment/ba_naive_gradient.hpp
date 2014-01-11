#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAICE_GRADIENT_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAICE_GRADIENT_HPP_

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
          int params_per_point>
struct BANaiveGradient
{
  typedef _Scalar Scalar;
  typedef _Index Index;
  typedef EIGEN_VECTOR(Scalar, params_per_camera) CameraSegment;
  typedef EIGEN_VECTOR(Scalar, params_per_point) PointSegment;
  typedef EIGEN_STD_VECTOR(CameraSegment) CameraSegmentContainer;
  typedef EIGEN_STD_VECTOR(PointSegment) PointSegmentContainer;

  Scalar operator[](Index i) const
  {
    Index camera_params_size = number_of_cameras * params_per_camera;
    Index params_size = camera_params_size +
                        number_of_points * params_per_point;
    if (i < 0 || i >= params_size)
    {
      return Scalar(0);
    }
    else if (i < camera_params_size)
    {
      return camera_segments[i / params_per_camera][i % params_per_camera];
    }
    else
    {
      Index off = i - camera_params_size;
      return point_segments[off / params_per_point][off % params_per_point];
    }
  }

  Index number_of_cameras;
  Index number_of_points;
  CameraSegmentContainer camera_segments;
  PointSegmentContainer point_segments;
};

}
}
}

#endif
