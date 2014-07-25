#ifndef _HS_SFM_CALIBRATE_PLANAR_CALIBRATOR_HPP_
#define _HS_SFM_CALIBRATE_PLANAR_CALIBRATOR_HPP_

#include "hs_sfm/sfm_utility/camera_type.hpp"

namespace hs
{
namespace sfm
{
namespace calibrate
{

template <typename _Scalar>
class PlanarCalibrator
{
public:
  typedef _Scalar Scalar;
  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef EIGEN_MATRIX(Scalar, 2, 2) KeyCovariance;
  typedef EIGEN_MATRIX(Scalar, 3, 3) PointCovariance;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef std::pair<Key, Point> Correspondence;
  typedef EIGEN_STD_VECTOR(Correspondence) PatternView;
  typedef EIGEN_STD_VECTOR(PatternView) PatternViewContainer;

public:
  int operator() (const PatternViewContainer& pattern_views,
                  const KeyCovariance& key_covariance,
                  const PointCovariance& point_covariance,
                  IntrinsicParams& intrinsic_params,
                  ExtrinsicParamsContainer& extrinsic_params_set) const
  {
    return 0;
  }
};

}
}
}

#endif
