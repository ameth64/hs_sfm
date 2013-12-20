#ifndef _HS_SFM_SFM_UTILITY_PROJECTIVE_FUNCTIONS_HPP_
#define _HS_SFM_SFM_UTILITY_PROJECTIVE_FUNCTIONS_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/radial_distortor.hpp"
#include "hs_sfm/sfm_utility/decentering_distortor.hpp"

namespace hs
{
namespace sfm
{

template <typename _Scalar>
class ProjectiveFunctions
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_VECTOR(Scalar, 2) Key;

  static Key WorldPointProjectToImageKey(
    const IntrinsicParams& intrinsic_params,
    const ExtrinsicParams& extrinsic_params,
    const Point3D& point)
  {
    Point3D camera_point = point - extrinsic_params.position();
    camera_point = extrinsic_params.rotation() * camera_point;
    camera_point /= camera_point[2];

    Key radial_delta_key;
    RadialDistortor<Scalar>()(intrinsic_params.k1(),
                              intrinsic_params.k2(),
                              intrinsic_params.k3(),
                              camera_point[0],
                              camera_point[1],
                              radial_delta_key[0],
                              radial_delta_key[1]);

    Key decentering_delta_key;
    DecenteringDistortor<Scalar>()(intrinsic_params.d1(),
                                   intrinsic_params.d2(),
                                   camera_point[0],
                                   camera_point[1],
                                   decentering_delta_key[0],
                                   decentering_delta_key[1]);

    return (camera_point.segment(0, 2) +
            radial_delta_key +
            decentering_delta_key);
  }

};

}
}

#endif
