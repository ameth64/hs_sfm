#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_ANALYTICAL_X_COVARIANCE_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_ANALYTICAL_X_COVARIANCE_CALCULATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_normal_equation_builder.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VectorFunction>
class CameraSharedAnalyticalXCovarianceCalculator;

/**
 *  解析式的计算多个相机共享的Bundle Adjustment向量函数的x向量协方差。
 *
 */
template <typename _Scalar>
class CameraSharedAnalyticalXCovarianceCalculator<
  CameraSharedVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;

  typedef CameraSharedAnalyticalJacobianMatrixCalculator<VectorFunction>
          JacobianMatrixCalculator;
  typedef typename JacobianMatrixCalculator::JacobianMatrix JacobianMatrix;

  typedef CameraSharedNormalEquationBuilder<Scalar> NormalEquationBuilder;
  typedef typename NormalEquationBuilder::NormalMatrix NormalMatrix;
  typedef typename NormalEquationBuilder::YCovarianceInverse
                   YCovarianceInverse;

  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic)
          DenseNormalMatrix;
  typedef DenseNormalMatrix XCovariance;

  Err operator() (const VectorFunction& vector_function,
                  const XVector& optimized_x,
                  const YCovarianceInverse& y_covariance_inverse,
                  XCovariance& x_covariance) const
  {
    JacobianMatrixCalculator jacobian;
    JacobianMatrix jacobian_matrix;
    if (jacobian(vector_function, optimized_x, jacobian_matrix) != 0)
      return -1;

    NormalEquationBuilder normal_equation_builder;
    NormalMatrix normal_matrix;
    if (normal_equation_builder.ComputeNormalMatrix(jacobian_matrix,
                                                    y_covariance_inverse,
                                                    normal_matrix) != 0)
      return -1;

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
    EIGEN_VECTOR(Scalar, Eigen::Dynamic) d = svd.singularValues();
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
