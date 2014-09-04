﻿#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_CERES_OPTIMIZOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_CERES_OPTIMIZOR_HPP_

#include <type_traits>

#include <ceres/problem.h>
#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

class CameraSharedProjection
{
public:
  CameraSharedProjection(
    double observed_x, double observed_y,
    double cov_inv_00, double cov_inv_01,
    double cov_inv_10, double cov_inv_11)
    : observed_x_(observed_x), observed_y_(observed_y),
      cov_inv_00_(cov_inv_00), cov_inv_01_(cov_inv_01),
      cov_inv_10_(cov_inv_10), cov_inv_11_(cov_inv_11)
  {
  }

  template <typename T>
  bool operator() (const T* const point,
                   const T* const image,
                   const T* const camera,
                   T* residuals) const
  {
    T camera_point[3];
    ceres::AngleAxisRotatePoint(image, point, camera_point);
    camera_point[0] += image[3];
    camera_point[1] += image[4];
    camera_point[2] += image[5];
    T normalized_x = camera_point[0] / camera_point[2];
    T normalized_y = camera_point[1] / camera_point[2];

    T r2 = normalized_x * normalized_x + normalized_y * normalized_y;
    T r4 = r2 * r2;
    T r6 = r2 * r4;
    T radial_coeff = camera[0] * r2 + camera[1] * r4 + camera[2] * r6;
    T radial_x = radial_coeff * normalized_x;
    T radial_y = radial_coeff * normalized_y;

    T decentering_x = 2.0 * camera[3] * normalized_x * normalized_y +
                      camera[4] * (r2 + 2.0 * normalized_x * normalized_x);
    T decentering_y = 2.0 * camera[4] * normalized_x * normalized_y +
                      camera[3] * (r2 + 2.0 * normalized_y * normalized_y);

    normalized_x += radial_x + decentering_x;
    normalized_y += radial_y + decentering_y;

    T image_x = camera[5] * normalized_x + camera[6] * normalized_y + camera[7];
    T image_y = camera[5] * camera[9] * normalized_y + camera[8];

    T diff_x = image_x - T(observed_x_);
    T diff_y = image_y - T(observed_y_);

    residuals[0] = T(cov_inv_00_) * diff_x + T(cov_inv_01_) * diff_y;
    residuals[1] = T(cov_inv_10_) * diff_x + T(cov_inv_11_) * diff_y;

    return true;
  }

  static ceres::CostFunction* Create(
    const double observed_x, const double observed_y,
    const double cov_inv_00, const double cov_inv_01,
    const double cov_inv_10, const double cov_inv_11)
  {
    return (new ceres::AutoDiffCostFunction<
                  CameraSharedProjection, 2, 3, 6, 10>(
              new CameraSharedProjection(observed_x, observed_y,
                                         cov_inv_00, cov_inv_01,
                                         cov_inv_10, cov_inv_11)));
  }
protected:
  double observed_x_;
  double observed_y_;

  double cov_inv_00_;
  double cov_inv_01_;
  double cov_inv_10_;
  double cov_inv_11_;
};

class CameraSharedConstraint
{
public:
  CameraSharedConstraint(double expected_value, double constraint)
    : expected_value_(expected_value), constraint_(constraint) {}

  template <typename T>
  bool operator() (const T* const value, T* residuals) const
  {
    residuals[0] = (value[0] - T(expected_value_)) * T(constraint_);
    return true;
  }

  static ceres::CostFunction* Create(double expected_value, double constraint)
  {
    return (new ceres::AutoDiffCostFunction<
                  CameraSharedConstraint, 1, 1>(
              new CameraSharedConstraint(expected_value, constraint)));
  }
private:
  double expected_value_;
  double constraint_;
};

template <typename _VectorFunction>
class CameraSharedCeresOptimizor;

template <typename _Scalar>
class CameraSharedCeresOptimizor<CameraSharedVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef CameraSharedYCovarianceInverse<Scalar> YCovarianceInverse;

  CameraSharedCeresOptimizor(const XVector& initial_x)
    : initial_x_(initial_x) {}

  Err operator() (const VectorFunction& vector_function,
                  const YVector& near_y,
                  const YCovarianceInverse& y_covariance_inverse,
                  XVector& optimized_x) const
  {
    const IntrinsicComputationsMask& intrinsic_computations_mask =
      vector_function.intrinsic_computations_mask();
    if (vector_function.is_fix_points() || vector_function.is_fix_images() ||
        vector_function.is_fix_cameras() ||
        !intrinsic_computations_mask[COMPUTE_RADIAL_DISTORTION] ||
        !intrinsic_computations_mask[COMPUTE_DECENTERING_DISTORTION] ||
        !intrinsic_computations_mask[COMPUTE_INTRINSIC_PARAMS])
    {
      return -1;
    }
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
    Index point_params_size = vector_function.GetPointParamsSize();
    Index extrinsic_params_size = vector_function.GetExtrinsicParamsSize();
    Index intrinsic_params_size_per_camera =
      vector_function.GetIntrinsicParamsSizePerCamera();
    Index x_point_begin = 0;
    Index x_image_begin = point_params_size;
    Index x_camera_begin = x_image_begin + extrinsic_params_size;
    for (Index i = 0; i < vector_function.number_of_keys(); i++)
    {
      Index image_id = vector_function.feature_maps()[i].first;
      Index point_id = vector_function.feature_maps()[i].second;
      Index camera_id = vector_function.image_camera_map()[image_id];

      double observed_x =
        double(near_y[i * VectorFunction::params_per_key_ + 0]);
      double observed_y =
        double(near_y[i * VectorFunction::params_per_key_ + 1]);

      const typename YCovarianceInverse::KeyBlock& key_block =
        y_covariance_inverse.GetKeyBlock(i);
      double cov_inv_00 = key_block(0, 0);
      double cov_inv_01 = key_block(0, 1);
      double cov_inv_10 = key_block(1, 0);
      double cov_inv_11 = key_block(1, 1);

      ceres::CostFunction* cost_function =
        CameraSharedProjection::Create(observed_x, observed_y,
                                       cov_inv_00, cov_inv_01,
                                       cov_inv_10, cov_inv_11);

      ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);

      double* point = x_data + x_point_begin +
                      point_id * VectorFunction::params_per_point_;
      double* image = x_data + x_image_begin +
                      image_id * VectorFunction::extrinsic_params_per_image_;
      double* camera = x_data + x_camera_begin +
                       camera_id * intrinsic_params_size_per_camera;

      problem.AddResidualBlock(
        cost_function, loss_function, point, image, camera);
    }

    Index y_offset = vector_function.GetYKeysSize();
    Index constraint_offset = 0;
    auto itr_point_constraint = vector_function.point_constraints().begin();
    auto itr_point_constraint_end = vector_function.point_constraints().end();
    for (; itr_point_constraint != itr_point_constraint_end;
         ++itr_point_constraint)
    {
      Index point_id = Index(itr_point_constraint->point_id);
      Index x_offset = point_id * VectorFunction::params_per_point_;
      if (itr_point_constraint->mask[POINT_CONSTRAIN_X])
      {
        double expected_value = double(near_y[y_offset]);
        double constraint =
          y_covariance_inverse.GetConstraint(constraint_offset);
        ceres::CostFunction* cost_function =
          CameraSharedConstraint::Create(expected_value, constraint);

        double* value = x_data + x_offset + 0;

        problem.AddResidualBlock(cost_function, NULL, value);

        y_offset++;
        constraint_offset++;
      }
      if (itr_point_constraint->mask[POINT_CONSTRAIN_Y])
      {
        double expected_value = double(near_y[y_offset]);
        double constraint =
          y_covariance_inverse.GetConstraint(constraint_offset);
        ceres::CostFunction* cost_function =
          CameraSharedConstraint::Create(expected_value, constraint);

        double* value = x_data + x_offset + 1;

        problem.AddResidualBlock(cost_function, NULL, value);

        y_offset++;
        constraint_offset++;
      }
      if (itr_point_constraint->mask[POINT_CONSTRAIN_Z])
      {
        double expected_value = double(near_y[y_offset]);
        double constraint =
          y_covariance_inverse.GetConstraint(constraint_offset);
        ceres::CostFunction* cost_function =
          CameraSharedConstraint::Create(expected_value, constraint);

        double* value = x_data + x_offset + 2;

        problem.AddResidualBlock(cost_function, NULL, value);

        y_offset++;
        constraint_offset++;
      }
    }

    auto itr_image_constraint = vector_function.image_constraints().begin();
    auto itr_image_constraint_end = vector_function.image_constraints().end();
    for (; itr_image_constraint != itr_image_constraint_end;
         ++itr_image_constraint)
    {
      Index image_id = Index(itr_image_constraint->image_id);
      Index x_offset = x_image_begin +
                       image_id * VectorFunction::extrinsic_params_per_image_;
      if (itr_image_constraint->mask[IMAGE_CONSTRAIN_ROTATION])
      {
        for (Index i = 0; i < 3; i++)
        {
          double expected_value = double(near_y[y_offset + i]);
          double constraint =
            y_covariance_inverse.GetConstraint(constraint_offset + i);
          ceres::CostFunction* cost_function =
            CameraSharedConstraint::Create(expected_value, constraint);

          double* value = x_data + x_offset + i;
          problem.AddResidualBlock(cost_function, NULL, value);
        }
        y_offset += 3;
        constraint_offset += 3;
      }
      for (Index i = 0; i < 3; i++)
      {
        if (itr_image_constraint->mask[IMAGE_CONSTRAIN_POSITION_X + i])
        {
          double expected_value = double(near_y[y_offset]);
          double constraint =
            y_covariance_inverse.GetConstraint(constraint_offset);
          ceres::CostFunction* cost_function =
            CameraSharedConstraint::Create(expected_value, constraint);

          double* value = x_data + x_offset + 3 + i;
          problem.AddResidualBlock(cost_function, NULL, value);

          y_offset++;
          constraint_offset++;
        }
      }
    }

    auto itr_camera_constraint = vector_function.camera_constraints().begin();
    auto itr_camera_constraint_end = vector_function.camera_constraints().end();
    for (; itr_camera_constraint != itr_camera_constraint_end;
         ++itr_camera_constraint)
    {
      const IntrinsicComputationsMask& intrinsic_computations_mask =
        vector_function.intrinsic_computations_mask();
      if ((itr_camera_constraint->radial_mask.any() &&
           !intrinsic_computations_mask[COMPUTE_RADIAL_DISTORTION]) ||
          (itr_camera_constraint->decentering_mask.any() &&
           !intrinsic_computations_mask[COMPUTE_DECENTERING_DISTORTION]) ||
          (itr_camera_constraint->intrinsic_mask.any() &&
           !intrinsic_computations_mask[COMPUTE_INTRINSIC_PARAMS]))
      {
        return -1;
      }
      Index camera_id = Index(itr_camera_constraint->camera_id);
      Index x_offset = x_camera_begin +
                       camera_id * intrinsic_params_size_per_camera;
      for (Index i = 0; i < 3; i++)
      {
        if (itr_camera_constraint->radial_mask[RADIAL_CONSTRAIN_K1 + i])
        {
          double expected_value = double(near_y[y_offset]);
          double constraint =
            y_covariance_inverse.GetConstraint(constraint_offset);
          ceres::CostFunction* cost_function =
            CameraSharedConstraint::Create(expected_value, constraint);

          double* value = x_data + x_offset + i;
          problem.AddResidualBlock(cost_function, NULL, value);

          y_offset++;
          constraint_offset++;
        }
      }
      if (intrinsic_computations_mask[COMPUTE_RADIAL_DISTORTION])
      {
        x_offset += 3;
      }

      for (Index i = 0; i < 2; i++)
      {
        if (itr_camera_constraint->radial_mask[DECENTERING_CONSTRAIN_D1 + i])
        {
          double expected_value = double(near_y[y_offset]);
          double constraint =
            y_covariance_inverse.GetConstraint(constraint_offset);
          ceres::CostFunction* cost_function =
            CameraSharedConstraint::Create(expected_value, constraint);

          double* value = x_data + x_offset + i;
          problem.AddResidualBlock(cost_function, NULL, value);

          y_offset++;
          constraint_offset++;
        }
      }
      if (intrinsic_computations_mask[COMPUTE_DECENTERING_DISTORTION])
      {
        x_offset += 2;
      }

      for (Index i = 0; i < 5; i++)
      {
        if (itr_camera_constraint->radial_mask[
              INTRINSIC_CONSTRAIN_FOCAL_LENGTH + i])
        {
          double expected_value = double(near_y[y_offset]);
          double constraint =
            y_covariance_inverse.GetConstraint(constraint_offset);
          ceres::CostFunction* cost_function =
            CameraSharedConstraint::Create(expected_value, constraint);

          double* value = x_data + x_offset + i;
          problem.AddResidualBlock(cost_function, NULL, value);

          y_offset++;
          constraint_offset++;
        }
      }
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout<<summary.FullReport()<<"\n";

    return 0;
  }

private:
  XVector initial_x_;
};

}
}
}

#endif
