#ifndef _HS_SFM_SFM_UTILITY_DEBUG_TMP_HPP_
#define _HS_SFM_SFM_UTILITY_DEBUG_TMP_HPP_

#include <string>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"

namespace hs
{
namespace sfm
{

template <typename _Scalar>
class DebugTrue
{
public:
  typedef _Scalar Scalar;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;

  DebugTrue(const IntrinsicParamsContainer& intrinsic_params_set_true =
              IntrinsicParamsContainer(),
            const ExtrinsicParamsContainer& extrinsic_params_set_true =
              ExtrinsicParamsContainer(),
            const Point3DContainer& points_true =
              Point3DContainer(),
            const std::string& prefix = std::string())
    : intrinsic_params_set_true_(intrinsic_params_set_true)
    , extrinsic_params_set_true_(extrinsic_params_set_true)
    , points_true_(points_true)
    , prefix_(prefix) {}

  const IntrinsicParamsContainer& intrinsic_params_set_true() const
  {
    return intrinsic_params_set_true_;
  }

  const ExtrinsicParamsContainer& extrinsic_params_set_true() const
  {
    return extrinsic_params_set_true_;
  }

  const Point3DContainer& points_true() const
  {
    return points_true_;
  }

  const std::string& prefix() const
  {
    return prefix_;
  }

private:
  const IntrinsicParamsContainer& intrinsic_params_set_true_;
  const ExtrinsicParamsContainer& extrinsic_params_set_true_;
  const Point3DContainer& points_true_;
  std::string prefix_;
};

}
}

#endif