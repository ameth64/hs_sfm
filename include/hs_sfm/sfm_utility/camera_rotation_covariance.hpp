#ifndef _HS_SFM_UTILITY_CAMERA_ROTATION_COVARIANCE_HPP_
#define _HS_SFM_UTILITY_CAMERA_ROTATION_COVARIANCE_HPP_

#include "hs_math/geometry/euler_angles.hpp"
#include "hs_math/geometry/rotation.hpp"
#include "hs_math/fdjac/forward_finite_difference_jacobian_matrix_calculator.hpp"
#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_math/linear_algebra/latraits/vector_eigen.hpp"
#include "hs_math/linear_algebra/latraits/matrix_eigen.hpp"
#include "hs_math/linear_algebra/lafunc/arithmetic_eigen.hpp"

namespace hs
{
namespace sfm
{

/**
 *  计算相机旋转的协方差。
 */
template <typename _Scalar>
class CameraRotaionCovarianceCalculator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::math::geometry::Rotation3D<Scalar> GeneralRotation;
  typedef hs::math::geometry::EulerAngles<Scalar> EulerAnglesRotation;

  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef typename Vector3::Index Index;

  struct GeneralToEulerAnglesVectorFunction
  {
    typedef EIGEN_VECTOR(Scalar, 3) XVector;
    typedef EIGEN_VECTOR(Scalar, 3) YVector;

    Err operator()(const XVector& x, YVector& y) const
    {
      Matrix33 ortho_matrix = GeneralRotation(x);
      EulerAnglesRotation euler_angles;
      euler_angles.template FromOrthoRotMat<2, 1, -3, 1>(ortho_matrix);
      y[0] = euler_angles[0];
      y[1] = euler_angles[1];
      y[2] = euler_angles[2];

      return 0;
    }
  };

  struct EulerAnglesToGeneralVectorFunction
  {
    typedef EIGEN_VECTOR(Scalar, 3) XVector;
    typedef EIGEN_VECTOR(Scalar, 3) YVector;

    Err operator()(const XVector& x, YVector& y) const
    {
      EulerAnglesRotation euler_angles(x[0], x[1], x[2]);
      Matrix33 ortho_matrix = 
        euler_angles.template ToOrthoRotMat<2, 1, -3, 1>();
      GeneralRotation general_rotation(ortho_matrix);
      y[0] = general_rotation[0];
      y[1] = general_rotation[1];
      y[2] = general_rotation[2];

      return 0;
    }
  };

  typedef hs::math::fdjac::ForwardFiniteDifferenceJacobianMatrixCalculator<
    EulerAnglesToGeneralVectorFunction> JacobianMatrixCalculator;
  typedef typename JacobianMatrixCalculator::JacobianMatrix JacobianMatrix;

  Err operator()(const Vector3& cameram_rotation,
                 Scalar x_rotation_stddev,
                 Scalar y_rotation_stddev,
                 Scalar z_rotation_stddev,
                 Matrix33& camera_rotation_covariance) const
  {
    GeneralToEulerAnglesVectorFunction general_to_euler_angles;
    Vector3 x;
    if (general_to_euler_angles(cameram_rotation, x) != 0)
    {
      return -1;
    }

    EulerAnglesToGeneralVectorFunction euler_angles_to_general;
    JacobianMatrixCalculator jacobian;
    JacobianMatrix jacobian_matrix;
    if (jacobian(euler_angles_to_general, x, jacobian_matrix) != 0)
    {
      return -1;
    }

    Matrix33 static_jacobian_matrix;
    for (Index i = 0; i < 3; i++)
    {
      for (Index j = 0; j < 3; j++)
      {
        static_jacobian_matrix(i, j) = jacobian_matrix(i, j);
      }
    }

    Matrix33 euler_angles_covariance = Matrix33::Zero();
    x_rotation_stddev *= Scalar(M_PI) / Scalar(180);
    y_rotation_stddev *= Scalar(M_PI) / Scalar(180);
    z_rotation_stddev *= Scalar(M_PI) / Scalar(180);
    euler_angles_covariance(0, 0) = x_rotation_stddev * x_rotation_stddev;
    euler_angles_covariance(1, 1) = y_rotation_stddev * y_rotation_stddev;
    euler_angles_covariance(2, 2) = z_rotation_stddev * z_rotation_stddev;

    camera_rotation_covariance = 
      (static_jacobian_matrix.transpose() *
       euler_angles_covariance.inverse() *
       static_jacobian_matrix).inverse();

    return 0;
  }

};

}
}

#endif
