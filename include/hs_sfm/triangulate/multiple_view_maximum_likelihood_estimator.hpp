#ifndef _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_MAXIMUM_LIKELIHOOD_ESTIMATOR_HPP_
#define _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_MAXIMUM_LIKELIHOOD_ESTIMATOR_HPP_

#include "hs_sfm/triangulate/multiple_view_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/triangulate/multiple_view_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace triangulate
{

template <typename _Scalar>
class MultipleViewMaximumLikelihoodEstimator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

private:
  typedef MultipleViewVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef MultipleViewLevenbergMarquardtOptimizor<VectorFunction> Optimizor;
  typedef typename Optimizor::XVector XVector;
  typedef typename Optimizor::YVector YVector;
  typedef typename Optimizor::YCovarianceInverse YCovarianceInverse;
  typedef typename hs::sfm::Undistortor<Scalar> Undistortor;

public:
  typedef typename VectorFunction::ExtrinsicParams ExtrinsicParams;
  typedef typename VectorFunction::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename VectorFunction::IntrinsicParams IntrinsicParams;
  typedef typename VectorFunction::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef EIGEN_STD_VECTOR(Key) KeyContainer;
  typedef XVector Point;

public:
  Err operator() (const IntrinsicParamsContainer& intrinsic_params_set,
                  const ExtrinsicParamsContainer& extrinsic_params_set,
                  const KeyContainer& keys,
                  Point& point) const
  {
    size_t number_of_cameras = intrinsic_params_set.size();
    if (extrinsic_params_set.size() != number_of_cameras ||
        keys.size() != number_of_cameras)
    {
      return -1;
    }
    if (number_of_cameras < 2)
    {
      return -1;
    }

    VectorFunction vector_function;
    vector_function.set_intrinsic_params_set(intrinsic_params_set);
    vector_function.set_extrinsic_params_set(extrinsic_params_set);

    Index y_size = Index(number_of_cameras * 2);
    YVector near_y(y_size);
    Undistortor undistortor;
    for (Index i = 0; i < Index(number_of_cameras); i++)
    {
      near_y.segment(i * VectorFunction::params_per_feature_,
                     VectorFunction::params_per_feature_) = keys[i];
    }

    YCovarianceInverse y_covariance_inverse(y_size, y_size);
    y_covariance_inverse.setIdentity();

    Optimizor optimizor;
    if (optimizor(vector_function, near_y, y_covariance_inverse, point) != 0)
    {
      return -1;
    }

    return 0;
  }
};

}
}
}

#endif
