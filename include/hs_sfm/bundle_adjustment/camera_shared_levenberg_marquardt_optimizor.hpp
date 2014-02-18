#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_LEVENBERG_MARQUARDT_OPTIMIZOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_LEVENBERG_MARQUARDT_OPTIMIZOR_HPP_

#include "hs_optimizor/nllso/levenberg_marquardt.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_nllso_meta_levenberg_marquardt.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VectorFunction>
class CameraSharedLevenbergMarquardtOptimizor;

template <typename _Scalar>
class CameraSharedLevenbergMarquardtOptimizor<
        CameraSharedVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef hs::optimizor::nllso::LevenbergMarquardt<VectorFunction> Optimizor;
  typedef typename Optimizor::XVector XVector;
  typedef typename Optimizor::YVector YVector;
  typedef typename Optimizor::YCovarianceInverse YCovarianceInverse;

  CameraSharedLevenbergMarquardtOptimizor(const XVector& initial_x)
    : initial_x_(initial_x), optimizor_() {}

  CameraSharedLevenbergMarquardtOptimizor(const XVector& initial_x,
                                          size_t max_itr_num,
                                          Scalar tau,
                                          Scalar eps1,
                                          Scalar eps2)
    : initial_x_(initial_x), optimizor_(max_itr_num, tau, eps1, eps2) {}

  Err operator() (const VectorFunction& vector_function,
                  const YVector& near_y,
                  const YCovarianceInverse& y_covariance_inverse,
                  XVector& optimized_x) const
  {
    Index x_size = vector_function.GetXSize();
    if (initial_x_.rows() != x_size)
    {
      return -1;
    }

    optimized_x.resize(x_size);
    optimized_x = initial_x_;
    return optimizor_(vector_function, near_y, y_covariance_inverse,
                      optimized_x);
  }

private:
  Optimizor optimizor_;
  XVector initial_x_;
};

}
}
}

#endif
