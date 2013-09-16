#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_GRADIENT_DELTA_X_DOT_PRODUCT_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_GRADIENT_DELTA_X_DOT_PRODUCT_CALCULATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_gradient.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BANaiveGradientDeltaXDotProductCalculator
{
public:
  typedef _Scalar Scalar;
  typedef BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector DeltaXVector;
  typedef BANaiveGradient<Scalar, Index,
                          VectorFunction::params_per_camera_,
                          VectorFunction::params_per_point_>
          Gradient;
  
  Scalar operator() (const Gradient& gradient,
                     const DeltaXVector& delta_x) const
  {
    Index rows = gradient.number_of_cameras *
                 VectorFunction::params_per_camera_ +
                 gradient.number_of_points *
                 VectorFunction::params_per_point_;

    Scalar dot_product = Scalar(0);
    for (Index i = 0; i < rows; i++)
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