#ifndef _HS_SFM_SFM_UTILITY_UNDISTORTOR_HPP_
#define _HS_SFM_SFM_UTILITY_UNDISTORTOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"

namespace hs
{
namespace sfm
{

template <typename _Scalar>
class Undistortor
{
public:
  typedef _Scalar Scalar;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_VECTOR(Scalar, 2) ImagePoint;

private:
  typedef typename IntrinsicParams::KMatrix KMatrix;
public:
  Undistortor(Scalar precision = Scalar(1e-8))
    : precision_(precision) {}

  Scalar precision() const {return precision_;}

  ImagePoint UndistortImagePoint(const ImagePoint& point_distort,
                                 const IntrinsicParams& intrinsic_params) const
  {
    ImagePoint point_undistort = point_distort;
    NormalizeImagePoint(intrinsic_params, point_undistort);
    UndistortNormalizedImagePoint(intrinsic_params, point_undistort);
    DenormalizedImagePoint(intrinsic_params, point_undistort);
    return point_undistort;
  }

private:
  void NormalizeImagePoint(const IntrinsicParams& intrinsic_params,
                           ImagePoint& point) const
  {
    KMatrix K_inv = intrinsic_params.GetKMatrix().inverse();
    point = K_inv.block(0, 0, 2, 2) * point + K_inv.block(0, 2, 2, 1);
  }

  void DenormalizedImagePoint(const IntrinsicParams& intrinsic_params,
                              ImagePoint& point) const
  {
    KMatrix K = intrinsic_params.GetKMatrix();
    point = K.block(0, 0, 2, 2) * point + K.block(0, 2, 2, 1);
  }

  void UndistortNormalizedImagePoint(const IntrinsicParams& intrinsic_params,
                                     ImagePoint& point) const
  {
    double k1 = double(intrinsic_params.k1());
    double k2 = double(intrinsic_params.k2());
    double k3 = double(intrinsic_params.k3());
    double d1 = double(intrinsic_params.d1());
    double d2 = double(intrinsic_params.d2());

    double x0 = double(point[0]);
    double y0 = double(point[1]);
    double x = double(point[0]);
    double y = double(point[1]);
    for (int i = 0; i < 5; i++)
    {
      double r2 = x * x + y * y;
      double radial_coff = 1.0 / (1 + r2 * (k1 + r2 * (k2 + r2 * k3)));
      double decentering_dx = 2.0 * d1 * x * y + d2 * (r2 + 2.0 * x * x);
      double decentering_dy = 2.0 * d2 * x * y + d1 * (r2 + 2.0 * y * y);
      x = (x0 - decentering_dx) * radial_coff;
      y = (y0 - decentering_dy) * radial_coff;
    }
    point[0] = Scalar(x);
    point[1] = Scalar(y);
  }

private:
  Scalar precision_;
};

}
}

#endif
