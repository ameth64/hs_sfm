#ifndef _HS_SFM_SFM_UTILITY_RADIAL_DISTORTOR_HPP_
#define _HS_SFM_SFM_UTILITY_RADIAL_DISTORTOR_HPP_

namespace hs
{
namespace sfm
{

template <typename _Scalar>
struct RadialDistortor
{
  typedef _Scalar Scalar;
  typedef int Err;
  Err operator() (Scalar k1, Scalar k2, Scalar k3,
                  Scalar x, Scalar y,
                  Scalar& dx, Scalar& dy) const
  {
    Scalar r2 = x * x + y * y;
    Scalar r4 = r2 * r2;
    Scalar r6 = r2 * r4;

    Scalar coeff = k1 * r2 + k2 * r4 + k3 * r6;

    dx = coeff * x;
    dy = coeff * y;

    return 0;
  }
};

}
}

#endif
