#ifndef _HS_SFM_SFM_UTILITY_DECENTERING_DISTORT_HPP_
#define _HS_SFM_SFM_UTILITY_DECENTERING_DISTORT_HPP_

namespace hs
{
namespace sfm
{

/**
 *  给定偏心畸变参数以及归一化相机座标系下的座标，计算偏心畸变偏移值。
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
 *        2d_1xy+d_2(r^2 + 2x^2) \\
 *        2d_2xy+d_1(r^2 + 2y^2)
 *      \end{array}
 *    \right]
 *  \f]
 *  其中\f$r^2 = x^2 + y^2\f$
 */
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
