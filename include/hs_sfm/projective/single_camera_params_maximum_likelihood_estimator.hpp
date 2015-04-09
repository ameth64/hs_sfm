#ifndef _HS_SFM_PROJECTIVE_SINGLE_CAMERA_PARAMS_MAXIMUM_LIKELIHOOD_ESTIMATOR_HPP_
#define _HS_SFM_PROJECTIVE_SINGLE_CAMERA_PARAMS_MAXIMUM_LIKELIHOOD_ESTIMATOR_HPP_

#include "hs_sfm/config/hs_config.hpp"

#if HS_HAVE_MKL
#define EIGEN_USE_MKL_ALL
#endif

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/projective/pmatrix_dlt_calculator.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_ceres_optimizor.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

template <typename _Scalar>
class SingleCameraParamsMaximumLikelihoodEstimator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_MATRIX(Scalar, 2, 2) KeyCovariance;

private:
  typedef PMatrixDLTCalculator<Scalar> DLTCalculator;
  typedef typename DLTCalculator::PMatrix PMatrix;
  typedef IntrinsicConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef IntrinsicConstrainedCeresOptimizor<VectorFunction>
          Optimizor;
  typedef typename Optimizor::XVector XVector;
  typedef typename Optimizor::YVector YVector;
  typedef typename Optimizor::YCovarianceInverse YCovarianceInverse;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;

  typedef hs::math::geometry::Rotation3D<Scalar> Rotation;

public:
  typedef typename DLTCalculator::Correspondence Correspondence;
  typedef typename DLTCalculator::CorrespondenceContainer
          CorrespondenceContainer;

public:
  Err operator() (const CorrespondenceContainer& correspondences,
                  const KeyCovariance& key_covariance,
                  const IntrinsicParams& intrinsic_params_initial,
                  Scalar focal_length_stddev,
                  Scalar skew_stddev,
                  Scalar principal_point_x_stddev,
                  Scalar principal_point_y_stddev,
                  Scalar pixel_ratio_stddev,
                  IntrinsicParams& intrinsic_params_estimate,
                  ExtrinsicParams& extrinsic_params_estimate) const
  {
    DLTCalculator dlt_calculator;
    PMatrix P;
    if (dlt_calculator(correspondences, P) != 0)
    {
      return -1;
    }

    Matrix33 K;
    Matrix33 R;
    Vector3 t;
    GetRTFromPMatrix(P, K, R, t);
    PMatrix P_test;
    P_test.template block<3, 3>(0, 0) = K * R;
    P_test.template block<3, 1>(0, 3) = K * t;

    YCovarianceInverse y_covariance_inverse;
    size_t number_of_correspondences = correspondences.size();
    GetYCovarianceInverse(key_covariance,
                          focal_length_stddev,
                          skew_stddev,
                          principal_point_x_stddev,
                          principal_point_y_stddev,
                          pixel_ratio_stddev,
                          number_of_correspondences,
                          y_covariance_inverse);

    VectorFunction vector_function;
    GetVectorFunction(correspondences, vector_function);

    YVector near_y;
    GetNearYVector(correspondences, intrinsic_params_initial, near_y);

    XVector initial_x;
    GetInitialXVector(intrinsic_params_initial, R, t, initial_x);

    //Optimizor optimizor(initial_x, 50, 1e-3, 1e-9, 1e-9);
    Optimizor optimizor(initial_x);
    XVector optimized_x;
    if (optimizor(vector_function, near_y, y_covariance_inverse,
                  optimized_x) != 0)
    {
      return -1;
    }

    GetIntrinsicAndExtrinsicParamsFromXVector(optimized_x,
                                              intrinsic_params_estimate,
                                              extrinsic_params_estimate);

    return 0;
  }

private:
  void GetRTFromPMatrix(const PMatrix& P,
                        Matrix33& K,
                        Matrix33& R,
                        Vector3& t) const
  {
    K = P.template block<3, 3>(0, 0);
    //对M作RQ分解
    Scalar nx = K.template block<1, 2>(2, 1).norm();
    Scalar cx = -K(2, 2) / nx;
    Scalar sx = K(2, 1) / nx;
    Matrix33 Qx;
    Qx << 1, 0, 0,
          0, cx, -sx,
          0, sx, cx;
    K = K * Qx;
    Scalar ny = std::sqrt(K(2, 0) * K(2, 0) + K(2, 2) * K(2, 2));
    Scalar cy = K(2, 2) / ny;
    Scalar sy = K(2, 0) / ny;
    Matrix33 Qy;
    Qy << cy, 0, sy,
          0, 1, 0,
          -sy, 0, cy;
    K = K * Qy;
    Scalar nz = K.template block<1, 2>(1, 0).norm();
    Scalar cz = -K(1, 1) / nz;
    Scalar sz = K(1, 0) / nz;
    Matrix33 Qz;
    Qz << cz, -sz, 0,
          sz, cz, 0,
          0, 0, 1;
    K = K * Qz;

    //使K矩阵对角线元素均为正
    PMatrix M = P;
    int negative = 0;
    if (K(0, 0) < Scalar(0)) negative++;
    if (K(1, 1) < Scalar(0)) negative++;
    if (K(2, 2) < Scalar(0)) negative++;
    if (negative % 2 == 1)
    {
      K *= -1;
      M *= -1;
    }
    Matrix33 fix = Matrix33::Identity();
    if (K(0, 0) < Scalar(0) && K(1, 1) < Scalar(0))
    {
      fix(0, 0) *= -1;
      fix(1, 1) *= -1;
    }
    else if (K(0, 0) < Scalar(0))
    {
      fix(0, 0) *= -1;
      fix(2, 2) *= -1;
    }
    else if (K(1, 1) < Scalar(0))
    {
      fix(1, 1) *= -1;
      fix(2, 2) *= -1;
    }

    K = K * fix;
    M = K.inverse() * M;
    R = M.template block<3, 3>(0, 0);
    t = M.col(3);

    K /= K(2, 2);
  }

  void GetYCovarianceInverse(const KeyCovariance& key_covariance,
                             Scalar focal_length_stddev,
                             Scalar skew_stddev,
                             Scalar principal_point_x_stddev,
                             Scalar principal_point_y_stddev,
                             Scalar pixel_ratio_stddev,
                             size_t number_of_correspondences,
                             YCovarianceInverse& y_covariance_inverse) const
  {
    y_covariance_inverse.key_blocks.resize(number_of_correspondences);
    KeyCovariance key_covariance_inverse = key_covariance.inverse();
    for (size_t i = 0; i < number_of_correspondences; i++)
    {
      y_covariance_inverse.key_blocks[i] = key_covariance_inverse;
    }
    y_covariance_inverse.focal_length_stddev =
      focal_length_stddev;
    y_covariance_inverse.skew_stddev =
      skew_stddev;
    y_covariance_inverse.principal_point_x_stddev =
      principal_point_x_stddev;
    y_covariance_inverse.principal_point_y_stddev =
      principal_point_y_stddev;
    y_covariance_inverse.pixel_ratio_stddev =
      pixel_ratio_stddev;
  }

  void GetVectorFunction(const CorrespondenceContainer& correspondences,
                         VectorFunction& vector_function) const
  {
    size_t number_of_points = correspondences.size();
    vector_function.points().resize(number_of_points);
    for (size_t i = 0; i < number_of_points; i++)
    {
      vector_function.points()[i] = correspondences[i].second;
    }
  }

  void GetNearYVector(const CorrespondenceContainer& correspondences,
                      const IntrinsicParams& intrinsic_params_initial,
                      YVector& near_y) const
  {
    Index number_of_keys = Index(correspondences.size());
    near_y.resize(number_of_keys * 2 + 5);
    for (Index i = 0; i < number_of_keys; i++)
    {
      near_y.segment(i * 2, 2) = correspondences[i].first;
    }
    near_y[number_of_keys * 2 + 0] =
      intrinsic_params_initial.focal_length();
    near_y[number_of_keys * 2 + 1] =
      intrinsic_params_initial.skew();
    near_y[number_of_keys * 2 + 2] =
      intrinsic_params_initial.principal_point_x();
    near_y[number_of_keys * 2 + 3] =
      intrinsic_params_initial.principal_point_y();
    near_y[number_of_keys * 2 + 4] =
      intrinsic_params_initial.pixel_ratio();
  }

  void GetInitialXVector(const Matrix33& K,
                         const Matrix33& R,
                         const Vector3& t,
                         XVector& initial_x) const
  {
    initial_x.resize(11);
    Rotation rotation = R;
    initial_x[0] = rotation[0];
    initial_x[1] = rotation[1];
    initial_x[2] = rotation[2];
    initial_x[3] = t[0];
    initial_x[4] = t[1];
    initial_x[5] = t[2];
    initial_x[6] = K(0, 0);
    initial_x[7] = K(0, 1);
    initial_x[8] = K(0, 2);
    initial_x[9] = K(1, 2);
    initial_x[10] = K(1, 1) / K(0, 0);
  }

  void GetInitialXVector(const IntrinsicParams& intrinsic_params_initial,
                         const Matrix33& R,
                         const Vector3& t,
                         XVector& initial_x) const
  {
    initial_x.resize(11);
    Rotation rotation = R;
    initial_x[0] = rotation[0];
    initial_x[1] = rotation[1];
    initial_x[2] = rotation[2];
    initial_x[3] = t[0];
    initial_x[4] = t[1];
    initial_x[5] = t[2];
    initial_x[6] = intrinsic_params_initial.focal_length();
    initial_x[7] = intrinsic_params_initial.skew();
    initial_x[8] = intrinsic_params_initial.principal_point_x();
    initial_x[9] = intrinsic_params_initial.principal_point_y();
    initial_x[10] = intrinsic_params_initial.pixel_ratio();
  }

  void GetIntrinsicAndExtrinsicParamsFromXVector(
    const XVector& x,
    IntrinsicParams& intrinsic_params_estimate,
    ExtrinsicParams& extrinsic_params_estimate) const
  {
    extrinsic_params_estimate.rotation()[0] = x[0];
    extrinsic_params_estimate.rotation()[1] = x[1];
    extrinsic_params_estimate.rotation()[2] = x[2];
    Vector3 t = x.template segment<3>(3);
    Matrix33 R = extrinsic_params_estimate.rotation();
    extrinsic_params_estimate.position() = -R.transpose() * t;

    intrinsic_params_estimate.set_focal_length(x[6]);
    intrinsic_params_estimate.set_skew(x[7]);
    intrinsic_params_estimate.set_principal_point_x(x[8]);
    intrinsic_params_estimate.set_principal_point_y(x[9]);
    intrinsic_params_estimate.set_pixel_ratio(x[10]);
  }

};

}
}
}

#endif
