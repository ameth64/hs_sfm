#ifndef _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_X_COVARIANCE_CALCULATOR_HPP_
#define _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_X_COVARIANCE_CALCULATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_math/fdjac/ffd_jac.hpp"

#include "hs_sfm/triangulate/multiple_view_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace triangulate
{

template <typename _VectorFunction>
class MultipleViewXCovarianceCalculator;

template <typename _Scalar>
class MultipleViewXCovarianceCalculator<MultipleViewVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef MultipleViewVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVec XVector;
  typedef typename VectorFunction::YVec YVector;

  typedef hs::math::fdjac::FwdFiniteDiffJacobian<VectorFunction> Jacobian;
  typedef typename Jacobian::Jac JacobianMatrix;

  typedef EIGEN_MAT(Scalar, Eigen::Dynamic, Eigen::Dynamic) YCovarianceInverse;
  typedef EIGEN_MAT(Scalar, VectorFunction::params_per_point_,
                            VectorFunction::params_per_point_) XCovariance;
  typedef typename XCovariance::Index Index;

  typedef int Err;

  Err operator() (const VectorFunction& vector_fucntion,
                  const XVector& optimized_x,
                  const YCovarianceInverse& y_covariance_inverse,
                  XCovariance& x_covariance) const
  {
    Jacobian jacobian;
    JacobianMatrix jacobian_matrix;
    if (jacobian(vector_fucntion, optimized_x, jacobian_matrix) != 0)
    {
      return -1;
    }

    XCovariance normal_matrix = jacobian_matrix.transpose() *
                                y_covariance_inverse *
                                jacobian_matrix;

    Eigen::JacobiSVD<XCovariance> svd(normal_matrix, 
                                      Eigen::ComputeFullU |
                                      Eigen::ComputeFullV);
    XCovariance D = XCovariance::Zero();
    EIGEN_VEC(Scalar, Eigen::Dynamic) d = svd.singularValues();
    const Scalar precision = Scalar(1e-6);
    for (Index i = 0; i < VectorFunction::params_per_point_; i++)
    {
      if (std::abs(d[i]) > precision)
      {
        D(i, i) = Scalar(1) / d[i];
      }
    }

    x_covariance = svd.matrixV() * D * (svd.matrixU().transpose());

    return 0;
  }
};

}
}
}

#endif