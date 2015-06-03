#ifndef _HS_SFM_PROJECTIVE_MLE_CERES_OPTIMIZOR_HPP_
#define _HS_SFM_PROJECTIVE_MLE_CERES_OPTIMIZOR_HPP_

#include <ceres/problem.h>
#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include "hs_sfm/projective/mle/intrinsic_constrained_vector_function.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_y_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

class IntrinsicConstrainedProjection
{
public:
  IntrinsicConstrainedProjection(
    const std::vector<double>& observes,
    const std::vector<double>& covariants,
    const std::vector<double>& points)
    : observes_(observes)
    , covariants_(covariants)
    , points_(points) {}

  template <typename T>
  bool operator() (const T* const x, T* residuals) const
  {
    size_t number_of_points = points_.size() / 3;
    if (number_of_points != observes_.size() / 2 ||
        number_of_points != covariants_.size() / 4) return false;

    for (size_t i = 0; i < number_of_points; i++)
    {
      T rotation[3];
      rotation[0] = x[0];
      rotation[1] = x[1];
      rotation[2] = x[2];
      T point[3];
      point[0] = T(points_[i * 3 + 0]);
      point[1] = T(points_[i * 3 + 1]);
      point[2] = T(points_[i * 3 + 2]);
      T camera_point[3];
      ceres::AngleAxisRotatePoint(rotation, point, camera_point);
      camera_point[0] += x[3];
      camera_point[1] += x[4];
      camera_point[2] += x[5];
      T normalized_x = camera_point[0] / camera_point[2];
      T normalized_y = camera_point[1] / camera_point[2];

      T r2 = normalized_x * normalized_x + normalized_y * normalized_y;
      T r4 = r2 * r2;
      T r6 = r2 * r4;
      T radial_coeff = x[11] * r2 + x[12] * r4 + x[13] * r6;
      T radial_x = radial_coeff * normalized_x;
      T radial_y = radial_coeff * normalized_y;

      T decentering_x = 2.0 * x[14] * normalized_x * normalized_y +
                        x[15] * (r2 + 2.0 * normalized_x * normalized_x);
      T decentering_y = 2.0 * x[15] * normalized_x * normalized_y +
                        x[14] * (r2 + 2.0 * normalized_y * normalized_y);

      T image_x = x[6] * normalized_x + x[7] * normalized_y + x[8];
      T image_y = x[6] * x[10] * normalized_y + x[9];
      T diff_x = image_x - T(observes_[i * 2 + 0]);
      T diff_y = image_y - T(observes_[i * 2 + 1]);

      T cov_inv_00 = T(covariants_[i * 4 + 0]);
      T cov_inv_01 = T(covariants_[i * 4 + 1]);
      T cov_inv_10 = T(covariants_[i * 4 + 2]);
      T cov_inv_11 = T(covariants_[i * 4 + 3]);

      residuals[i * 2 + 0] = T(cov_inv_00) * diff_x + T(cov_inv_01) * diff_y;
      residuals[i * 2 + 1] = T(cov_inv_10) * diff_x + T(cov_inv_11) * diff_y;

    }

    return true;
  }

  static ceres::CostFunction* Create(
    const std::vector<double>& observes,
    const std::vector<double>& covariants,
    const std::vector<double>& points)
  {
    return (new ceres::AutoDiffCostFunction<
                  IntrinsicConstrainedProjection, ceres::DYNAMIC, 16>(
              new IntrinsicConstrainedProjection(observes, covariants, points),
              int(observes.size())));
  }

protected:
  const std::vector<double>& observes_;
  const std::vector<double>& covariants_;
  const std::vector<double>& points_;
};

template <typename _VectorFunction>
class IntrinsicConstrainedCeresOptimizor;

template <typename _Scalar>
class IntrinsicConstrainedCeresOptimizor<
        IntrinsicConstrainedVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef IntrinsicConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef IntrinsicConstrainedYCovarianceInverse<Scalar> YCovarianceInverse;

private:
  typedef typename VectorFunction::Index Index;
  typedef typename YCovarianceInverse::KeyBlock KeyBlock;

public:
  IntrinsicConstrainedCeresOptimizor(const XVector& initial_x)
    : initial_x_(initial_x) {}

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

    optimized_x = initial_x_;
    double* x_data = nullptr;
    EIGEN_VECTOR(double, Eigen::Dynamic) casted_x;
    if (!std::is_same<Scalar, double>::value)
    {
      casted_x = optimized_x.template cast<double>();
      x_data = casted_x.data();
    }
    else
    {
      x_data = optimized_x.data();
    }

    ceres::Problem problem;
    size_t number_of_points = vector_function.points().size();
    double expected_focal_length = near_y[Index(number_of_points) * 2 + 0];
    double expected_skew = near_y[Index(number_of_points) * 2 + 1];
    double expected_principal_x = near_y[Index(number_of_points) * 2 + 2];
    double expected_principal_y = near_y[Index(number_of_points) * 2 + 3];
    double expected_pixel_ratio = near_y[Index(number_of_points) * 2 + 4];
    double expected_k1 = near_y[Index(number_of_points) * 2 + 5];
    double expected_k2 = near_y[Index(number_of_points) * 2 + 6];
    double expected_k3 = near_y[Index(number_of_points) * 2 + 7];
    double expected_d1 = near_y[Index(number_of_points) * 2 + 8];
    double expected_d2 = near_y[Index(number_of_points) * 2 + 9];

    std::vector<double> observes;
    std::vector<double> covariants;
    std::vector<double> points;
    for (size_t i = 0; i < number_of_points; i++)
    {
      observes.push_back(double(near_y[i * 2 + 0]));
      observes.push_back(double(near_y[i * 2 + 1]));
      const KeyBlock& key_block = y_covariance_inverse.key_blocks[i];
      covariants.push_back(double(key_block(0, 0)));
      covariants.push_back(double(key_block(0, 1)));
      covariants.push_back(double(key_block(1, 0)));
      covariants.push_back(double(key_block(1, 1)));
      points.push_back(double(vector_function.points()[i][0]));
      points.push_back(double(vector_function.points()[i][1]));
      points.push_back(double(vector_function.points()[i][2]));
    }
    ceres::CostFunction* cost_function =
      IntrinsicConstrainedProjection::Create(observes, covariants, points);
    problem.AddResidualBlock(cost_function, NULL, x_data);
    problem.SetParameterLowerBound(
      x_data, 6, expected_focal_length -
                 y_covariance_inverse.focal_length_stddev);
    problem.SetParameterUpperBound(
      x_data, 6, expected_focal_length +
                 y_covariance_inverse.focal_length_stddev);
    problem.SetParameterLowerBound(
      x_data, 7, expected_skew -
                 y_covariance_inverse.skew_stddev);
    problem.SetParameterUpperBound(
      x_data, 7, expected_skew +
                 y_covariance_inverse.skew_stddev);
    problem.SetParameterLowerBound(
      x_data, 8, expected_principal_x -
                 y_covariance_inverse.principal_point_x_stddev);
    problem.SetParameterUpperBound(
      x_data, 8, expected_principal_x +
                 y_covariance_inverse.principal_point_x_stddev);
    problem.SetParameterLowerBound(
      x_data, 9, expected_principal_y -
                 y_covariance_inverse.principal_point_y_stddev);
    problem.SetParameterUpperBound(
      x_data, 9, expected_principal_y +
                 y_covariance_inverse.principal_point_y_stddev);
    problem.SetParameterLowerBound(
      x_data, 10, expected_pixel_ratio -
                 y_covariance_inverse.pixel_ratio_stddev);
    problem.SetParameterUpperBound(
      x_data, 10, expected_pixel_ratio +
                 y_covariance_inverse.pixel_ratio_stddev);
    problem.SetParameterLowerBound(
      x_data, 11, expected_k1 -
                  y_covariance_inverse.k1_stddev);
    problem.SetParameterUpperBound(
      x_data, 11, expected_k1 +
                  y_covariance_inverse.k1_stddev);
    problem.SetParameterLowerBound(
      x_data, 12, expected_k2 -
                  y_covariance_inverse.k2_stddev);
    problem.SetParameterUpperBound(
      x_data, 12, expected_k2 +
                  y_covariance_inverse.k2_stddev);
    problem.SetParameterLowerBound(
      x_data, 13, expected_k3 -
                  y_covariance_inverse.k3_stddev);
    problem.SetParameterUpperBound(
      x_data, 13, expected_k3 +
                  y_covariance_inverse.k3_stddev);
    problem.SetParameterLowerBound(
      x_data, 14, expected_d1 -
                  y_covariance_inverse.d1_stddev);
    problem.SetParameterUpperBound(
      x_data, 14, expected_d1 +
                  y_covariance_inverse.d1_stddev);
    problem.SetParameterLowerBound(
      x_data, 15, expected_d2 -
                  y_covariance_inverse.d2_stddev);
    problem.SetParameterUpperBound(
      x_data, 15, expected_d2 +
                  y_covariance_inverse.d2_stddev);

    ceres::Solver::Options options;
    options.max_num_iterations = 50;
    options.function_tolerance = 1e-9;
    options.parameter_tolerance = 1e-9;
    options.logging_type = ceres::SILENT;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    if (!std::is_same<Scalar, double>::value)
    {
      optimized_x = casted_x.template cast<Scalar>();
    }

    return 0;
  }

private:
  XVector initial_x_;
};

}
}
}

#endif
