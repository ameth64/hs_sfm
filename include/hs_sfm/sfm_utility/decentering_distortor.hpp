#ifndef _HS_SFM_SFM_UTILITY_DECENTERING_DISTORT_HPP_
#define _HS_SFM_SFM_UTILITY_DECENTERING_DISTORT_HPP_

namespace hs
{
namespace sfm
{

template <typename _Scalar>
struct DecenteringDistortor
{
  typedef _Scalar Scalar;
  typedef int Err;
  Err operator() (Scalar d1, Scalar d2,
                  Scalar x, Scalar y,
                  Scalar& dx, Scalar& dy) const
  {
    Scalar r2 = x * x + y * y;

    dx = Scalar(2) * d1 * x * y + d2 * (r2 + Scalar(2) * x * x);
    dy = Scalar(2) * d2 * x * y + d1 * (r2 + Scalar(2) * y * y);

    return -1;
  }
};

}
}

#endif
