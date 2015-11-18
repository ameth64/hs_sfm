#ifndef _HS_SFM_PROJECTIVE_SINGLE_CAMERA_EXTRINSIC_PARAMS_OPTIMIZOR_HPP_
#define _HS_SFM_PROJECTIVE_SINGLE_CAMERA_EXTRINSIC_PARAMS_OPTIMIZOR_HPP_

#include <ceres/problem.h>
#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include "hs_sfm/config/hs_config.hpp"

#if HS_HAVE_MKL
#define EIGEN_USE_MKL_ALL
#endif

#include "hs_sfm/sfm_utility/camera_type.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

class SingleCameraProjection
{
public:
  SingleCameraProjection(const std::vector<double>& observes,
                         const std::vector<double>& points,
                         const std::vector<double>& intrinsic_params)
    : observes_(observes)
    , points_(points)
    , intrinsic_params_(intrinsic_params) {}

  template <typename T>
  bool operator() (const T* const x, T* residuals) const
  {
    size_t number_of_points = points_.size() / 3;
    if (number_of_points != observes_.size() / 2 ||
        intrinsic_params_.size() != 10) return false;

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

      T focal_length = T(intrinsic_params_[0]);
      T skew = T(intrinsic_params_[1]);
      T principal_x = T(intrinsic_params_[2]);
      T principal_y = T(intrinsic_params_[3]);
      T pixel_ratio = T(intrinsic_params_[4]);
      T k1 = T(intrinsic_params_[5]);
      T k2 = T(intrinsic_params_[6]);
      T k3 = T(intrinsic_params_[7]);
      T d1 = T(intrinsic_params_[8]);
      T d2 = T(intrinsic_params_[9]);

      T r2 = normalized_x * normalized_x + normalized_y * normalized_y;
      T r4 = r2 * r2;
      T r6 = r2 * r4;
      T radial_coeff = k1 * r2 + k2 * r4 + k3 * r6;
      T radial_x = radial_coeff * normalized_x;
      T radial_y = radial_coeff * normalized_y;

      T decentering_x = 2.0 * d1 * normalized_x * normalized_y +
                        d2 * (r2 + 2.0 * normalized_x * normalized_x);
      T decentering_y = 2.0 * d2 * normalized_x * normalized_y +
                        d1 * (r2 + 2.0 * normalized_y * normalized_y);

      T image_x = focal_length * normalized_x + skew * normalized_y +
                  principal_x;
      T image_y = focal_length * pixel_ratio * normalized_y + principal_y;
      T diff_x = image_x - T(observes_[i * 2 + 0]);
      T diff_y = image_y - T(observes_[i * 2 + 1]);

      residuals[i * 2 + 0] = diff_x;
      residuals[i * 2 + 1] = diff_y;
    }

    return true;
  }

  static ceres::CostFunction* Create(
    const std::vector<double>& observes,
    const std::vector<double>& points,
    const std::vector<double>& intrinsic_params)
  {
    return (new ceres::AutoDiffCostFunction<
                  SingleCameraProjection, ceres::DYNAMIC, 6>(
              new SingleCameraProjection(observes, points, intrinsic_params),
              int(observes.size())));
  }

private:
  const std::vector<double>& observes_;
  const std::vector<double>& points_;
  const std::vector<double>& intrinsic_params_;
};

/**
 *  给定三维点和二维点的对应，内参数，以及初始的外参数，计算优化后的外参数。
 *
 */
template <typename _Scalar>
class SingleCameraExtrinsicParamsOptimizor
{
public:
  typedef _Scalar Scalar;
  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef std::pair<Key, Point> Correspondence;
  typedef EIGEN_STD_VECTOR(Correspondence) CorrespondenceContainer;

private:
  typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;

public:
  SingleCameraExtrinsicParamsOptimizor() {}

public:
  int operator() (const CorrespondenceContainer& correspondences,
                  const IntrinsicParams& intrinsic_params,
                  ExtrinsicParams& extrinsic_params) const
  {
    size_t number_of_points = correspondences.size();
    std::vector<double> observes(number_of_points * 2);
    std::vector<double> points(number_of_points * 3);
    std::vector<double> intrinsic_params_optimize(10);
    double x[6];
    Point t = -(extrinsic_params.rotation() * extrinsic_params.position());
    x[0] = double(extrinsic_params.rotation()[0]);
    x[1] = double(extrinsic_params.rotation()[1]);
    x[2] = double(extrinsic_params.rotation()[2]);
    x[3] = double(t[0]);
    x[4] = double(t[1]);
    x[5] = double(t[2]);
    for (size_t i = 0; i < number_of_points; i++)
    {
      observes[i * 2 + 0] = double(correspondences[i].first[0]);
      observes[i * 2 + 1] = double(correspondences[i].first[1]);
      points[i * 3 + 0] = double(correspondences[i].second[0]);
      points[i * 3 + 1] = double(correspondences[i].second[1]);
      points[i * 3 + 2] = double(correspondences[i].second[2]);
    }
    intrinsic_params_optimize[0] = double(intrinsic_params.focal_length());
    intrinsic_params_optimize[1] = double(intrinsic_params.skew());
    intrinsic_params_optimize[2] = double(intrinsic_params.principal_point_x());
    intrinsic_params_optimize[3] = double(intrinsic_params.principal_point_y());
    intrinsic_params_optimize[4] = double(intrinsic_params.pixel_ratio());
    intrinsic_params_optimize[5] = double(intrinsic_params.k1());
    intrinsic_params_optimize[6] = double(intrinsic_params.k2());
    intrinsic_params_optimize[7] = double(intrinsic_params.k3());
    intrinsic_params_optimize[8] = double(intrinsic_params.d1());
    intrinsic_params_optimize[9] = double(intrinsic_params.d2());

    ceres::CostFunction* cost_function =
      SingleCameraProjection::Create(
        observes, points, intrinsic_params_optimize);
    ceres::Problem problem;
    problem.AddResidualBlock(cost_function, NULL, x);
    ceres::Solver::Options options;
    options.max_num_iterations = 50;
    options.function_tolerance = 1e-9;
    options.parameter_tolerance = 1e-9;
    options.logging_type = ceres::SILENT;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    extrinsic_params.rotation()[0] = x[0];
    extrinsic_params.rotation()[1] = x[1];
    extrinsic_params.rotation()[2] = x[2];
    RMatrix rmatrix = extrinsic_params.rotation();
    t << Scalar(x[3]),
         Scalar(x[4]),
         Scalar(x[5]);
    extrinsic_params.position() = -rmatrix.transpose() * t;

    return 0;
  }

};

}
}
}

#endif
