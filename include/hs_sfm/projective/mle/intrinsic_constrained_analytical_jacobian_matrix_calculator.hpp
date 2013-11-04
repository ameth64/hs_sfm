#ifndef _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_ANALYTICAL_JACOBIAN_MATRIX_CALCULATOR_HPP_
#define _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_ANALYTICAL_JACOBIAN_MATRIX_CALCULATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/projective/mle/intrinsic_constrained_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

template <typename _VectorFunction>
class IntrinsicConstrainedAnalyticalJacobianMatrixCalculator;

template <typename _Scalar>
class IntrinsicConstrainedAnalyticalJacobianMatrixCalculator<
        IntrinsicConstrainedVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  
  typedef IntrinsicConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;

  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic) JacobianMatrix;

private:
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

public:
  Err operator() (const VectorFunction& vector_function,
                  const XVector& x,
                  JacobianMatrix& jacobian_matrix) const
  {
    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    jacobian_matrix.resize(y_size, x_size);
    jacobian_matrix.setZero();
    Vector3 r = x.template segment<3>(0);
    Vector3 t = x.template segment<3>(3);
    Scalar theta = r.norm();
    Scalar theta2 = theta * theta;
    Scalar theta4 = theta2 * theta2;
    Matrix33 K;
    K << x[6], x[7], x[8],
         0, x[6] * x[10], x[9],
         0, 0, 1;

    Index number_of_points = Index(vector_function.points().size());
    for (Index i = 0; i < number_of_points; i++)
    {
      Vector3 p = vector_function.points()[i];
      Scalar d = r.dot(p);

      //sin(theta)
      Scalar st = sin(theta);
      //cos(theta)
      Scalar ct = cos(theta);

      //f = cos(theta) * p
      Vector3 f = ct * p;
      //\frac{\partial f}{\partial r}
      Matrix33 pfpr = p * (-st / theta * r).transpose();

      //s = sin(theta) / theta
      Scalar s = st / theta;
      //\frac{\partial s}{\partial r}
      Vector3 pspr = (ct * theta - st) / theta2 / theta * r;

      //v = r x p
      Vector3 v = r.cross(p);
      //\frac{\partial v}{\partial r}
      Matrix33 pvpr;
      pvpr << 0, p[2], -p[1],
              -p[2], 0, p[0],
              p[1], -p[0], 0;

      //c = (1 - cos(theta)) / theta^2
      Scalar c = (1 - ct) / theta2;
      //\frac{partial c}{\partial r}
      Vector3 pcpr =
        (st * theta - 2 * (1 - ct)) / theta4 * r;

      //u = r^T * p * r
      Vector3 u = d * r;
      //\frac{partial u}{\partial r}
      Matrix33 pupr;
      pupr << d + r[0] * p[0], r[0] * p[1], r[0] * p[2],
              r[1] * p[0], d + r[1] * p[1], r[1] * p[2],
              r[2] * p[0], r[2] * p[1], d + r[2] * p[2];

      Vector3 g = f + s * v + c * u + t;
      //\frac{\partial g}{\partial r}
      Matrix33 pgpr = pfpr +
                      v * pspr.transpose() + s * pvpr +
                      u * pcpr.transpose() + c * pupr;
      Matrix33 pgpt = Matrix33::Identity();
      Vector3 y = K * g;
      Matrix33 pypr = K * pgpr;
      Matrix33 pypt = K * pgpt;

      for (Index m = 0; m < 3; m++)
      {
        jacobian_matrix(i * 2, m) =
          (pypr(0, m) * y[2] - y[0] * pypr(2, m)) / y[2] / y[2];
        jacobian_matrix(i * 2 + 1, m) =
          (pypr(1, m) * y[2] - y[1] * pypr(2, m)) / y[2] / y[2];

        jacobian_matrix(i * 2, 3 + m) =
          (pypt(0, m) * y[2] - y[0] * pypt(2, m)) / y[2] / y[2];
        jacobian_matrix(i * 2 + 1, 3 + m) =
          (pypt(1, m) * y[2] - y[1] * pypt(2, m)) / y[2] / y[2];
      }

      jacobian_matrix(i * 2, 6) = g[0] / y[2];
      jacobian_matrix(i * 2 + 1, 6) = g[1] * x[10] / y[2];
      jacobian_matrix(i * 2, 7) = g[1] / y[2];
      jacobian_matrix(i * 2 + 1, 7) = 0;
      jacobian_matrix(i * 2, 8) = g[2] / y[2];
      jacobian_matrix(i * 2 + 1, 8) = 0;
      jacobian_matrix(i * 2, 9) = 0;
      jacobian_matrix(i * 2 + 1, 9) = g[2] / y[2];
      jacobian_matrix(i * 2, 10) = 0;
      jacobian_matrix(i * 2 + 1, 10) = g[1] * x[6] / y[2];
    }

    jacobian_matrix.template block<5, 5>(number_of_points * 2, 6).setIdentity();

    return 0;
  }
};

}
}
}

#endif