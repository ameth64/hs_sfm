#ifndef _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_REDISUALS_HPP_
#define _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_REDISUALS_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/triangulate/multiple_view_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace triangulate
{

template <typename _VectorFunction>
class MultipleViewResidualsCalculator;

template <typename _Scalar>
class MultipleViewResidualsCalculator<MultipleViewVectorFunction<_Scalar> >
{

public:
  typedef _Scalar Scalar;
  typedef MultipleViewVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVec XVector;
  typedef typename VectorFunction::YVec YVector;
  typedef EIGEN_MAT(Scalar, Eigen::Dynamic, Eigen::Dynamic) YCovInv;
  typedef EIGEN_VEC(Scalar, Eigen::Dynamic) Residuals;
  
  typedef int Err;
  typedef typename Residuals::Index Index;

  Err operator() (const VectorFunction& vector_function,
                  const XVector& noised_x,
                  const YVector& estimated_y,
                  Residuals& r) const
  {
    YVector noised_y;
    if (vector_function(noised_x, noised_y)) return -1;

    if (estimated_y.rows() != noised_y.rows()) return -1;

    r = noised_y - estimated_y;

    return 0;
  }

  Residuals operator()(const YVector& noised_y, const YVector& estimated_y) const
  {
    return noised_y - estimated_y;
  }

  Scalar operator()(const Residuals& residuals, 
                    const YCovInv& y_covariance_inverse) const
  {
    return residuals.transpose() * y_covariance_inverse * residuals;
  }
};

}
}
}

#endif