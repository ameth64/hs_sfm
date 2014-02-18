#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_GRADIENT_DELTA_X_DOT_PRODUCT_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_GRADIENT_DELTA_X_DOT_PRODUCT_CALCULATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_gradient.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedGradientDeltaXDotProductCalculator
{
public:
  typedef _Scalar Scalar;
  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector DeltaXVector;
  typedef CameraSharedGradient<Scalar> Gradient;

  Scalar operator() (const Gradient& gradient,
                     const DeltaXVector& delta_x) const
  {
    Index x_size = gradient.GetXSize();
    Scalar dot_product = Scalar(0);
    for (Index i = 0; i < x_size; i++)
    {
      dot_product += gradient[i] * delta_x[i];
    }
    return dot_product;
  }
};

}
}
}

#endif
