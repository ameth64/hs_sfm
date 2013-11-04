#ifndef _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_NORMAL_MATRIX_BUILDER_HPP_
#define _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_NORMAL_MATRIX_BUILDER_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/projective/mle/intrinsic_constrained_y_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

template <typename _Scalar>
class IntrinsicConstrainedNormalEquationBuilder
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic) JacobianMatrix;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) Residuals;
  typedef IntrinsicConstrainedYCovarianceInverse<Scalar> YCovarianceInverse;
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic) NormalMatrix;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) Gradient;

private:
  typedef typename NormalMatrix::Index Index;

public:
  Err operator() (const JacobianMatrix& jacobian_matrix,
                  const Residuals& residuals,
                  const YCovarianceInverse& y_covariance_inverse,
                  NormalMatrix& normal_matrix,
                  Gradient& gradient) const
  {
    Index number_of_points = Index(y_covariance_inverse.key_blocks.size());
    Index x_size = jacobian_matrix.cols();
    normal_matrix.resize(x_size, x_size);
    gradient.resize(x_size);
    for (Index i = 0; i < x_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        Scalar key_sum = Scalar(0);
        for (Index k = 0; k < number_of_points; k++)
        {
          key_sum += (jacobian_matrix.block(k * 2, i, 2, 1).transpose() *
                      y_covariance_inverse.key_blocks[k] *
                      jacobian_matrix.block(k * 2, j, 2, 1))(0, 0);
        }

        Scalar intrinsic_sum =
          jacobian_matrix(number_of_points * 2 + 0, i) *
          jacobian_matrix(number_of_points * 2 + 0, j) /
          (y_covariance_inverse.focal_length_stddev *
           y_covariance_inverse.focal_length_stddev) +
          jacobian_matrix(number_of_points * 2 + 1, i) *
          jacobian_matrix(number_of_points * 2 + 1, j) /
          (y_covariance_inverse.skew_stddev *
           y_covariance_inverse.skew_stddev) +
          jacobian_matrix(number_of_points * 2 + 2, i) *
          jacobian_matrix(number_of_points * 2 + 2, j) /
          (y_covariance_inverse.principal_point_x_stddev *
           y_covariance_inverse.principal_point_x_stddev) +
          jacobian_matrix(number_of_points * 2 + 3, i) *
          jacobian_matrix(number_of_points * 2 + 3, j) /
          (y_covariance_inverse.principal_point_y_stddev *
           y_covariance_inverse.principal_point_y_stddev) +
          jacobian_matrix(number_of_points * 2 + 4, i) *
          jacobian_matrix(number_of_points * 2 + 4, j) /
          (y_covariance_inverse.pixel_ratio_stddev *
           y_covariance_inverse.pixel_ratio_stddev);

        normal_matrix(i, j) = key_sum + intrinsic_sum;
      }
    }

    for (Index i = 0; i < x_size; i++)
    {
      Scalar key_sum = Scalar(0);
      for (Index j = 0; j < number_of_points; j++)
      {
        key_sum += (jacobian_matrix.block(j * 2, i, 2, 1).transpose() *
                    y_covariance_inverse.key_blocks[j] *
                    residuals.segment(j * 2, 2))[0];
      }

      Scalar intrinsic_sum =
        jacobian_matrix(number_of_points * 2 + 0, i) *
        residuals[number_of_points * 2 + 0] /
        (y_covariance_inverse.focal_length_stddev *
          y_covariance_inverse.focal_length_stddev) +
        jacobian_matrix(number_of_points * 2 + 1, i) *
        residuals[number_of_points * 2 + 1] /
        (y_covariance_inverse.skew_stddev *
          y_covariance_inverse.skew_stddev) +
        jacobian_matrix(number_of_points * 2 + 2, i) *
        residuals[number_of_points * 2 + 2] /
        (y_covariance_inverse.principal_point_x_stddev *
          y_covariance_inverse.principal_point_x_stddev) +
        jacobian_matrix(number_of_points * 2 + 3, i) *
        residuals[number_of_points * 2 + 3] /
        (y_covariance_inverse.principal_point_y_stddev *
          y_covariance_inverse.principal_point_y_stddev) +
        jacobian_matrix(number_of_points * 2 + 4, i) *
        residuals[number_of_points * 2 + 4] /
        (y_covariance_inverse.pixel_ratio_stddev *
          y_covariance_inverse.pixel_ratio_stddev);

      gradient[i] = key_sum + intrinsic_sum;
    }

    return 0;
  }
};

}
}
}

#endif