#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NOISED_X_GENERATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NOISED_X_GENERATOR_HPP_

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/utility/camera_rotation_covariance.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BANaiveNoisedXGenerator
{
public:
  typedef int Err;
  typedef _Scalar Scalar;
  typedef BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;

private:
  typedef typename VectorFunction::Index Index;
  typedef hs::sfm::CameraRotaionCovarianceCalculator<Scalar>
    CameraRotationCovariance;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;

public:
  Err operator()(const VectorFunction& vector_function,
                 const XVector& true_x,
                 Scalar x_rotation_stddev,
                 Scalar y_rotation_stddev,
                 Scalar z_rotation_stddev,
                 Scalar x_pos_stddev,
                 Scalar y_pos_stddev,
                 Scalar z_pos_stddev,
                 Scalar x_point_stddev,
                 Scalar y_point_stddev,
                 Scalar z_point_stddev,
                 XVector& noised_x) const
  {
    noised_x.resize(true_x.rows());
    noised_x = true_x;
    //put some noise to near_x
    Index number_of_camera = vector_function.number_of_cameras();
    Index camera_params_size = vector_function.GetCameraParamsSize();
    Index x_params_size = vector_function.GetXSize();
    CameraRotationCovariance rotation_covariance_calculator;
    for (Index i = 0; i < number_of_camera; i++)
    {
      Vector3 camera_rotation = true_x.segment(
        i * VectorFunction::params_per_camera_, 3);
      Matrix33 camera_rotation_covariance;
      if (rotation_covariance_calculator(camera_rotation, 
        x_rotation_stddev,
        y_rotation_stddev,
        z_rotation_stddev,
        camera_rotation_covariance) != 0)
      {
        std::cout<<"rotation covariance calculator failed!\n";
        return -1;
      }
      Vector3 noised_camera_rotation;
      hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
        camera_rotation, camera_rotation_covariance, noised_camera_rotation);
      noised_x.segment(i * VectorFunction::params_per_camera_, 3) =
        noised_camera_rotation;

      Vector3 camera_pos = true_x.segment(
        i * VectorFunction::params_per_camera_ + 3, 3);
      Matrix33 camera_pos_covariance = Matrix33::Zero();
      camera_pos_covariance(0, 0) = x_pos_stddev * x_pos_stddev;
      camera_pos_covariance(1, 1) = y_pos_stddev * y_pos_stddev;
      camera_pos_covariance(2, 2) = z_pos_stddev * z_pos_stddev;
      Vector3 noised_camera_pos;
      hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
        camera_pos, camera_pos_covariance, noised_camera_pos);
      noised_x.segment(i * VectorFunction::params_per_camera_ + 3, 3) =
        noised_camera_pos;
    }

    Index number_of_point = vector_function.number_of_points();
    for (Index i = 0; i < number_of_point; i++)
    {
      Vector3 point = true_x.segment(camera_params_size + 
        i * VectorFunction::params_per_point_, 3);
      Matrix33 point_covariance = Matrix33::Zero();
      point_covariance(0, 0) = x_point_stddev * x_point_stddev;
      point_covariance(1, 1) = y_point_stddev * y_point_stddev;
      point_covariance(2, 2) = z_point_stddev * z_point_stddev;
      Vector3 noised_point;
      hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
        point, point_covariance, noised_point);
      noised_x.segment(camera_params_size + 
        i * VectorFunction::params_per_point_, 3) = noised_point;
    }

    return 0;
  }
};

}
}
}

#endif
