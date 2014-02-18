#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_GRADIENT_INFINITE_NORM_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_GRADIENT_INFINITE_NORM_CALCULATOR_HPP_

#include <algorithm>
#include <limits>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_gradient.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedGradientInfiniteNormCalculator
{
public:
  typedef _Scalar Scalar;
  typedef CameraSharedGradient<Scalar> Gradient;
  typedef typename Gradient::Index Index;

  Scalar operator() (const Gradient& gradient) const
  {
    Index x_size = gradient.GetXSize();
    Scalar infinite_norm = -std::numeric_limits<Scalar>::max();
    for (Index i = 0; i < x_size; i++)
    {
      infinite_norm = std::max(gradient[i], infinite_norm);
    }

    return infinite_norm;
  }
};

}
}
}

#endif
