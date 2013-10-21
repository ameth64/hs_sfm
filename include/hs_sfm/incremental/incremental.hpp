#ifndef _HS_SFM_INCREMENTAL_INCREMENTAL_HPP_
#define _HS_SFM_INCREMENTAL_INCREMENTAL_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar>
class IncrementalSFM
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef ImageKeys<Scalar> Keys;
  typedef EIGEN_STD_VECTOR(Keys) KeysContainer;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;

public:
  Err operator()() const
  {
    return 0;
  }

private:
  Err RefineMatches() const
  {
    return 0;
  }

  Err GenerateTracks() const
  {
    return 0;
  }

  Err CalibrateInitialPair() const
  {
    return 0;
  }

  Err IncrementalAddCamera() const
  {
    return 0;
  }

private:
  KeysContainer keys_set_;
  MatchContainer refined_matches_;
  TrackContainer tracks_;
  size_t identity_camera_id_;
  size_t relative_camera_id_;
  ExtrinsicParamsContainer extrinsic_params_set_;
  Point3DContainer points_;
};

}
}
}

#endif