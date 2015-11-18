#ifndef _HS_SFM_SFM_UTILITY_RADIAL_DISTORTOR_HPP_
#define _HS_SFM_SFM_UTILITY_RADIAL_DISTORTOR_HPP_

namespace hs
{
namespace sfm
{

/**
 *  给定径向畸变参数以及归一化相机座标系下的座标，计算径向畸变偏移值。
 *
 *  \f[
 *    \left[
 *      \begin{array}{l}
 *        dx \\
 *        dy
 *      \end{array}
 *    \right] =
 *    \left[
 *      \begin{array}{l}
 *        (k_1r^2+k_2r^4+k_3r^6)x \\
 *        (k_1r^2+k_2r^4+k_3r^6)y \\
 *      \end{array}
 *    \right]
 *  \f]
 *  其中\f$r^2 = x^2 + y^2\f$
 */
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
