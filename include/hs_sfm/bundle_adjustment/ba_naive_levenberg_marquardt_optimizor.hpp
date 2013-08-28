#ifndef _HS_SFM_BUNDLE_ADJUSTEMENT_BA_NAIVE_LEVENBERG_MARQUARDT_OPTIMIZOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTEMENT_BA_NAIVE_LEVENBERG_MARQUARDT_OPTIMIZOR_HPP_

#include "hs_optimizor/nllso/levenberg_marquardt.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_nllso_meta.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VectorFunction>
class BANaiveLevenbergMarquardtOptimizor;

template <typename _Scalar>
class BANaiveLevenbergMarquardtOptimizor<BANaiveVecFunc<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef BANaiveVecFunc<_Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef hs::optimizor::nllso::LevenbergMarquardt<VectorFunction> Optimizor;
  typedef typename Optimizor::XVec XVector;
  typedef typename Optimizor::YVec YVector;
  typedef typename Optimizor::YCovInv YCovarianceInverse;
  typedef typename Optimizor::YCovInv YCovInv;

  BANaiveLevenbergMarquardtOptimizor(const XVector& initial_x)
    : initial_x_(initial_x), optimizor_() {}

  BANaiveLevenbergMarquardtOptimizor(const XVector& initial_x,
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
    Index x_size = vector_function.getXSize();
    if (initial_x_.rows() != x_size)
    {
      return -1;
    }

    optimized_x.resize(x_size);
    optimized_x = initial_x_;
    Err rst = optimizor_(vector_function, near_y, y_covariance_inverse,
                         optimized_x);
#ifdef _DEBUG
    YVector optimized_y;
    vector_function(optimized_x, optimized_y);
    YVector optimized_residuals = optimized_y - near_y;
    Scalar optimized_norm = optimized_residuals.squaredNorm();
    optimized_norm /= optimized_residuals.rows();
    optimized_norm = std::sqrt(optimized_norm);
    optimized_norm *= 2951.1378875492117335208598481088;
    
    YVector initial_y;
    vector_function(initial_x_, initial_y);
    YVector initial_residuals = initial_y - near_y;
    Scalar initial_norm = initial_residuals.squaredNorm();
    initial_norm /= initial_residuals.rows();
    initial_norm = std::sqrt(initial_norm);
    initial_norm *= 2951.1378875492117335208598481088;
#endif

    return rst;
  }

private:
  Optimizor optimizor_;
  XVector initial_x_;
};

}
}
}

#endif
