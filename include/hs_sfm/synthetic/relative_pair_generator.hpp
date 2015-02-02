#ifndef _HS_SFM_SYNTHETIC_RELATIVE_PAIR_GENERATOR_HPP_
#define _HS_SFM_SYNTHETIC_RELATIVE_PAIR_GENERATOR_HPP_

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"
#include "hs_math/geometry/euler_angles.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{
namespace synthetic
{

template <typename _Scalar, typename _ImageDimension>
class RelativePairGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef typename ExtrinsicParams::Position Position;
  typedef typename ExtrinsicParams::Rotation Rotation;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef ImageParams<ImageDimension> Image;
  typedef EIGEN_STD_VECTOR(Image) ImageContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;

private:
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

public:
  RelativePairGenerator(
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    const IntrinsicParams& intrinsic_params_identity,
    const IntrinsicParams& intrinsic_params_relative,
    Scalar lateral_overlap_ratio,
    Scalar longitudinal_overlap_ratio,
    Scalar scene_max_height,
    Scalar camera_height_stddev,
    Scalar camera_planar_stddev,
    Scalar camera_rotation_stddev)
    : ground_resolution_(ground_resolution),
      image_width_(image_width),
      image_height_(image_height),
      pixel_size_(pixel_size),
      intrinsic_params_identity_(intrinsic_params_identity),
      intrinsic_params_relative_(intrinsic_params_relative),
      lateral_overlap_ratio_(lateral_overlap_ratio),
      longitudinal_overlap_ratio_(longitudinal_overlap_ratio),
      scene_max_height_(scene_max_height),
      camera_height_stddev_(camera_height_stddev),
      camera_planar_stddev_(camera_planar_stddev),
      camera_rotation_stddev_(camera_rotation_stddev) {}

  Err GenerateExtrinsicParams(ExtrinsicParams& extrinsic_params_identity,
                              ExtrinsicParams& extrinsic_params_relative,
                              Image& image_identity,
                              Image& image_relative) const
  {
    image_identity.m_width = image_width_;
    image_identity.m_height = image_height_;
    image_identity.m_id = 0;
    image_identity.m_path = "synthetic";

    image_relative.m_width = image_width_;
    image_relative.m_height = image_height_;
    image_relative.m_id = 1;
    image_identity.m_path = "synthetic";

    Matrix33 camera_position_covariance = Matrix33::Identity();
    camera_position_covariance(0, 0) = camera_planar_stddev_;
    camera_position_covariance(1, 1) = camera_planar_stddev_;
    camera_position_covariance(2, 2) = camera_height_stddev_;
    Matrix33 camera_rotation_covariance = Matrix33::Identity();
    camera_rotation_covariance *= camera_rotation_stddev_;

    Position mean_position_identity = Position::Zero();
    mean_position_identity(2) = intrinsic_params_identity_.focal_length() *
                                ground_resolution_;
    hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
      mean_position_identity,
      camera_position_covariance,
      extrinsic_params_identity.position());
    Vector3 mean_angles = Vector3::Zero();
    Vector3 angles_identity;
    hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
      mean_angles, camera_rotation_covariance, angles_identity);
    hs::math::geometry::EulerAngles<Scalar> euler_angles_identity(
      angles_identity[0] / 180 * Scalar(M_PI),
      angles_identity[1] / 180 * Scalar(M_PI),
      angles_identity[2] / 180 * Scalar(M_PI));
    Matrix33 R_identity =
      euler_angles_identity.template ToOrthoRotMat<2, 1, -3, 1>();
    R_identity.transposeInPlace();
    R_identity.col(2) *= -1;
    R_identity.col(1) *= -1;
    extrinsic_params_identity.rotation() = R_identity;

    Position mean_position_relative;
    mean_position_relative <<
      (Scalar(1) - lateral_overlap_ratio_) * Scalar(image_width_) *
      ground_resolution_,
      (Scalar(1) - longitudinal_overlap_ratio_) * Scalar(image_height_) *
      ground_resolution_,
      intrinsic_params_relative_.focal_length() * ground_resolution_;
    hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
      mean_position_relative,
      camera_position_covariance,
      extrinsic_params_relative.position());

    Vector3 angles_relative;
    hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
      mean_angles, camera_rotation_covariance, angles_relative);
    hs::math::geometry::EulerAngles<Scalar> euler_angles_relative(
      angles_relative[0] / 180 * Scalar(M_PI),
      angles_relative[1] / 180 * Scalar(M_PI),
      angles_relative[2] / 180 * Scalar(M_PI));
    Matrix33 R_relative =
      euler_angles_relative.template ToOrthoRotMat<2, 1, -3, 1>();
    R_relative.transposeInPlace();
    R_relative.col(2) *= -1;
    R_relative.col(1) *= -1;
    extrinsic_params_relative.rotation() = R_relative;

    return 0;
  }

  Err GeneratePoints(size_t number_of_points,
                     Point3DContainer& points) const
  {
    points.resize(number_of_points);
    Point3D min_point;
    min_point <<
      (Scalar(0.5) - lateral_overlap_ratio_) * image_width_ *
      ground_resolution_,
      (Scalar(0.5) - longitudinal_overlap_ratio_) * image_height_ *
      ground_resolution_,
      0;
    Point3D max_point;
    max_point <<
      Scalar(0.5) * Scalar(image_width_) * ground_resolution_,
      Scalar(0.5) * Scalar(image_height_) * ground_resolution_,
      scene_max_height_;

    for (size_t i = 0; i < number_of_points; i++)
    {
      hs::math::random::UniformRandomVar<Scalar, 3>::Generate(
        min_point, max_point, points[i]);
    }

    return 0;
  }

  const IntrinsicParams& intrinsic_params_identity() const
  {
    return intrinsic_params_identity_;
  }

  const IntrinsicParams& intrinsic_params_relative() const
  {
    return intrinsic_params_relative_;
  }

  Scalar ground_resolution() const
  {
    return ground_resolution_;
  }

  ImageDimension image_width() const
  {
    return image_width_;
  }

  ImageDimension image_height() const
  {
    return image_height_;
  }

private:
  Scalar ground_resolution_;
  ImageDimension image_width_;
  ImageDimension image_height_;
  Scalar pixel_size_;
  IntrinsicParams intrinsic_params_identity_;
  IntrinsicParams intrinsic_params_relative_;
  Scalar lateral_overlap_ratio_;
  Scalar longitudinal_overlap_ratio_;
  Scalar scene_max_height_;
  Scalar camera_height_stddev_;
  Scalar camera_planar_stddev_;
  Scalar camera_rotation_stddev_;
};

}
}
}

#endif
