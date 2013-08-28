#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_XCOV_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_XCOV_CALCULATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_analytic_jac.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_nllso_normal_equation_builder.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VectorFunction>
class BANaiveXCovCalculator;

template <typename _Scalar>
class BANaiveXCovCalculator<BANaiveVecFunc<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef BANaiveVecFunc<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVec XVector;
  typedef typename VectorFunction::YVec YVecctor;

  typedef BANaiveAnalyticJacobian<VectorFunction> Jacobian;
  typedef typename Jacobian::Jac JacobianMatrix;

  typedef BANaiveNormalEquationBuilder<Scalar> NormalEquationBuilder;
  typedef typename NormalEquationBuilder::NormalMat NormalMatrix;
  typedef typename NormalEquationBuilder::YCovInv YCovarianceInverse;

  typedef EIGEN_MAT(Scalar, Eigen::Dynamic, Eigen::Dynamic) DenseNormalMatrix;
  typedef DenseNormalMatrix XCovariance;

  typedef int Err;

  Err operator()(const VectorFunction& vector_function,
                 const XVector& optimized_x,
                 const YCovarianceInverse& y_covariance_inverse,
                 XCovariance& x_covariance) const
  {
    Jacobian jacobian;
    JacobianMatrix jacobian_matrix;
    if (jacobian(vector_function, optimized_x, jacobian_matrix) != 0)
    {
      return -1;
    }

    NormalEquationBuilder normal_equation_builder;
    NormalMatrix normal_matrix;
    if (normal_equation_builder.BuildNormalMatrix(jacobian_matrix,
                                                  y_covariance_inverse,
                                                  normal_matrix) != 0)
    {
      return -1;
    }

    Index x_size = optimized_x.rows();

    DenseNormalMatrix dense_normal_matrix(x_size, x_size);
    for (Index i = 0; i < x_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        dense_normal_matrix(i, j) = normal_matrix.coeff(i, j);
      }
    }

    //计算正规矩阵的伪逆即为X的协方差矩阵
    Eigen::JacobiSVD<DenseNormalMatrix> svd(dense_normal_matrix,
                                            Eigen::ComputeThinU |
                                            Eigen::ComputeThinV);
    DenseNormalMatrix D = DenseNormalMatrix::Zero(x_size, x_size);
    EIGEN_VEC(Scalar, Eigen::Dynamic) d = svd.singularValues();
    const Scalar precision = Scalar(1e-6);
    for (Index i = 0; i < x_size; i++)
    {
      if (std::abs(d(i)) > precision)
      {
        D(i, i) = Scalar(1) / d(i);
      }
    }
    x_covariance.resize(x_size, x_size);
    x_covariance = svd.matrixV() * D * (svd.matrixU().transpose());

    return 0;
  }
};

}
}
}

#endif
