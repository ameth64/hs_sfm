#ifndef _HS_SFM_SFM_PIPELINE_GLOBAL_TYPES_HPP_
#define _HS_SFM_SFM_PIPELINE_GLOBAL_TYPES_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_math/geometry/rotation.hpp"

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
struct RotationPair
{
  typedef _Scalar Scalar;
  typedef size_t Index;
  typedef hs::math::geometry::Rotation3D<Scalar> Rotation;

  Index first_id;
  Index second_id;
  Rotation rotation;
};

template <typename _Scalar>
struct PositionPair
{
  typedef _Scalar Scalar;
  typedef size_t Index;
  typedef EIGEN_VECTOR(Scalar, 3) Position;

  Index first_id;
  Index second_id;
  Position position;
  Scalar wieght;
};

}
}
}

#endif
