#ifndef _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_LEVERNBERG_MARQUARDT_OPTIMIZOR_HPP_
#define _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_LEVERNBERG_MARQUARDT_OPTIMIZOR_HPP_

#include "hs_optimizor/nllso/meta_eigen_levenberg_marquardt.hpp"
#include "hs_optimizor/nllso/levenberg_marquardt.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/undistortor.hpp"
#include "hs_sfm/triangulate/multiple_view_vector_function.hpp"
#include "hs_sfm/triangulate/multiple_view_dlt.hpp"

namespace hs
{
namespace sfm
{
namespace triangulate
{

template <typename _VectorFunction>
class MultipleViewLevenbergMarquardtOptimizor;

template <typename _Scalar>
class MultipleViewLevenbergMarquardtOptimizor<
        MultipleViewVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef MultipleViewVectorFunction<Scalar> VectorFunction;
  typedef hs::optimizor::nllso::LevenbergMarquardt<VectorFunction> Optimizor;
  typedef typename Optimizor::XVector XVector;
  typedef typename Optimizor::YVector YVector;
  typedef typename Optimizor::YCovarianceInverse YCovarianceInverse;

private:
  typedef typename VectorFunction::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename VectorFunction::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename XVector::Index Index;

  typedef MultipleViewDLT<Scalar> Initializor;
  typedef typename Initializor::PMatrix PMatrix;
  typedef typename Initializor::PMatrixContainer PMatrixContainer;
  typedef typename Initializor::HomogeneousKey HomogeneousKey;
  typedef typename Initializor::HomogeneousKeyContainer
          HomogeneousKeyContainer;
  typedef hs::sfm::CameraFunctions<Scalar> Camera;
  typedef typename VectorFunction::ExtrinsicParams ExtrinsicParams;
  typedef typename ExtrinsicParams::Rotation Rotation;
  typedef typename ExtrinsicParams::Position Position;

public:
  Err operator()(const VectorFunction& vector_function,
                 const YVector& near_y,
                 const YCovarianceInverse& y_covariance_inverse,
                 XVector& optimized_x) const
  {
    typedef EIGEN_VECTOR(Scalar, 2) Key;
    typedef hs::sfm::Undistortor<Scalar> Undistortor;
    typedef typename Undistortor::IntrinsicParams IntrinsicParams;

    const IntrinsicParamsContainer& intrins =
      vector_function.intrinsic_params_set();
    const ExtrinsicParamsContainer& extrins =
      vector_function.extrinsic_params_set();
    Index number_of_views = intrins.size();
    if (number_of_views != extrins.size() || number_of_views < 2)
    {
      return -1;
    }

    PMatrixContainer pmatrices;
    HomogeneousKeyContainer homogeneous_keys;
    Undistortor undistortor;
    for (Index i = 0; i < number_of_views; i++)
    {
      PMatrix pmatrix = Camera::GetProjectionMatrix(intrins[i], extrins[i]);
      pmatrices.push_back(pmatrix);
      HomogeneousKey homogeneous_key;
      Key key_distort =
        near_y.segment(i * VectorFunction::params_per_feature_,
                       VectorFunction::params_per_feature_);
      Key key_undistort =
        undistortor.UndistortImagePoint(key_distort, intrins[i]);
      homogeneous_key.segment(0, 2) = key_undistort;
      homogeneous_key[2] = Scalar(1);
      homogeneous_keys.push_back(homogeneous_key);
    }

    Initializor initializor;
    if (initializor(pmatrices, homogeneous_keys, optimized_x) != 0) return -1;

    Optimizor optimizor;
    if (optimizor(vector_function, near_y, y_covariance_inverse,
                  optimized_x) != 0) return -1;

    return 0;
  }
};

}
}
}

#endif
