#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_GRADIENT_INFINITE_NORM_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_GRADIENT_INFINITE_NORM_CALCULATOR_HPP_

#include <algorithm>
#include <limits>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_gradient.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BANaiveGradientInfiniteNormCalculator
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

  Scalar operator() (const Gradient& gradient) const
  {
    Index rows = gradient.number_of_cameras *
                 VectorFunction::params_per_camera_ +
                 gradient.number_of_points *
                 VectorFunction::params_per_point_;

    Scalar infinite_norm = -std::numeric_limits<Scalar>::max();
    
    for (Index i = 0; i < rows; i++)
    {
      infinite_norm = std::max(std::abs(gradient[i]), infinite_norm);
    }

    return infinite_norm;
  }
};

}
}
}

#endif
