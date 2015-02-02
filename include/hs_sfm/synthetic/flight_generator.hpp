#ifndef _HS_SFM_SYNTHETIC_FLIGHT_GENERATOR_HPP_
#define _HS_SFM_SYNTHETIC_FLIGHT_GENERATOR_HPP_

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
class FlightGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef typename ExtrinsicParams::Position Position;
  typedef typename ExtrinsicParams::Rotation Rotation;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef ImageParams<ImageDimension> Image;
  typedef EIGEN_STD_VECTOR(Image) ImageContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;

private:
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

public:
  FlightGenerator(Scalar focal_length_in_metre,
                  size_t number_of_strips,
                  size_t number_of_cameras_in_strip,
                  Scalar ground_resolution,
                  ImageDimension image_width,
                  ImageDimension image_height,
                  Scalar pixel_size,
                  size_t number_of_points,
                  Scalar lateral_overlap_ratio,
                  Scalar longitudinal_overlap_ratio,
                  Scalar scene_max_height,
                  Scalar camera_height_stddev,
                  Scalar camera_planar_stddev,
                  Scalar camera_rotation_stddev)
       : focal_length_in_metre_(focal_length_in_metre),
         number_of_strips_(number_of_strips),
         number_of_cameras_in_strip_(number_of_cameras_in_strip), 
         ground_resolution_(ground_resolution),
         image_width_(image_width), image_height_(image_height),
         pixel_size_(pixel_size),
         number_of_points_(number_of_points),
         lateral_overlap_ratio_(lateral_overlap_ratio),
         longitudinal_overlap_ratio_(longitudinal_overlap_ratio),
         scene_max_height_(scene_max_height),
         camera_height_stddev_(camera_height_stddev),
         camera_planar_stddev_(camera_planar_stddev),
         camera_rotation_stddev_(camera_rotation_stddev) {}

  Scalar GetFocalLengthInPixel() const
  {
    return focal_length_in_metre_ / pixel_size_;
  }

  Err operator ()(ExtrinsicParamsContainer& extrinsic_params_set,
                  ImageContainer& images,
                  Point3DContainer& points) const
  {
    GenerateCameras(extrinsic_params_set, images);
    GeneratePoints(number_of_points_, points);

    return 0;
  }

  inline Scalar GetFlightHeight() const
  {
    return focal_length_in_metre_ * ground_resolution_ / pixel_size_;
  }

  void GetDimensions(Scalar& scene_x_dimension,
                     Scalar& scene_y_dimension,
                     Scalar& scene_z_dimension) const
  {
    Scalar longitudinal_range_per_camera =
      (1 - longitudinal_overlap_ratio_) * ground_resolution_ *
      Scalar(image_height_);
    Scalar lateral_range_per_camera = 
      (1 - lateral_overlap_ratio_) * ground_resolution_ *
      Scalar(image_width_);

    scene_x_dimension = image_width_ * ground_resolution_ + 
      lateral_range_per_camera * (number_of_strips_ - 1);
    scene_y_dimension = image_height_ * ground_resolution_ +
      longitudinal_range_per_camera * (number_of_cameras_in_strip_ - 1);
    scene_z_dimension = scene_max_height_;
  }

  void GeneratePoints(size_t number_of_points,
                      Point3DContainer& points) const
  {
    points.resize(number_of_points);
    Scalar scene_x_dimension, scene_y_dimension, scene_z_dimension;
    GetDimensions(scene_x_dimension, scene_y_dimension, scene_z_dimension);

    Point3D max_point;
    max_point << scene_x_dimension * 0.5,
                 scene_y_dimension * 0.5,
                 scene_z_dimension;
    Point3D min_point;
    min_point << -scene_x_dimension * 0.5,
                 -scene_y_dimension * 0.5,
                 0;

    for (size_t i = 0; i < number_of_points; i++)
    {
      hs::math::random::UniformRandomVar<Scalar, 3>::Generate(
        min_point, max_point, points[i]);
    }
  }

  void GeneratePlanarPoints(size_t number_of_points,
                            Scalar height,
                            Point3DContainer& points) const
  {
    points.resize(number_of_points);
    Scalar scene_x_dimension, scene_y_dimension, scene_z_dimension;
    GetDimensions(scene_x_dimension, scene_y_dimension, scene_z_dimension);

    EIGEN_VECTOR(Scalar, 2) max_planar_point;
    max_planar_point << scene_x_dimension * 0.5,
                        scene_y_dimension * 0.5;
    EIGEN_VECTOR(Scalar, 2) min_planar_point = -max_planar_point;
    for (size_t i = 0; i < number_of_points; i++)
    {
      EIGEN_VECTOR(Scalar, 2) planar_point;
      hs::math::random::UniformRandomVar<Scalar, 2>::Generate(
        min_planar_point, max_planar_point, planar_point);
      points[i].template segment<2>(0) = planar_point;
      points[i][2] = height;
    }
  }

  void GenerateCameras(ExtrinsicParamsContainer& extrinsic_params_set,
                       ImageContainer& images) const
  {
    Scalar flight_height = GetFlightHeight();

    Scalar longitudinal_range_per_camera =
      (1 - longitudinal_overlap_ratio_) * ground_resolution_ *
      Scalar(image_height_);
    Scalar lateral_range_per_camera = 
      (1 - lateral_overlap_ratio_) * ground_resolution_ *
      Scalar(image_width_);

    size_t number_of_cameras = number_of_strips_ * number_of_cameras_in_strip_;

    extrinsic_params_set.resize(number_of_cameras);
    images.resize(number_of_cameras);

    Scalar scene_x_dimension, scene_y_dimension, scene_z_dimension;
    GetDimensions(scene_x_dimension, scene_y_dimension, scene_z_dimension);

    Matrix33 camera_position_covariance = Matrix33::Identity();
    camera_position_covariance(0, 0) = camera_planar_stddev_;
    camera_position_covariance(1, 1) = camera_planar_stddev_;
    camera_position_covariance(2, 2) = camera_height_stddev_;
    Matrix33 camera_rotation_covariance = Matrix33::Identity();
    camera_rotation_covariance *= camera_rotation_stddev_;
    for (size_t i = 0; i < number_of_strips_; i++)
    {
      for (size_t j = 0; j < number_of_cameras_in_strip_; j++)
      {
        size_t id = i * number_of_cameras_in_strip_ + j;
        images[id].m_width = image_width_;
        images[id].m_height = image_height_;
        images[id].m_id = id;
        images[id].m_path = "test";

        Position mean_camera_position;
        mean_camera_position << -scene_x_dimension * 0.5 +
                                image_width_ * ground_resolution_ * 0.5 +
                                i * lateral_range_per_camera,
                                -scene_y_dimension * 0.5 +
                                image_height_ * ground_resolution_ * 0.5 + 
                                j * longitudinal_range_per_camera,
                                flight_height;

        hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
          mean_camera_position,
          camera_position_covariance,
          extrinsic_params_set[id].position());
        Position meanAngles = Position::Zero();
        Position angles;
        hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
          meanAngles, camera_rotation_covariance, angles);

        hs::math::geometry::EulerAngles<Scalar> euler_angles(
          angles[0] / 180 * Scalar(M_PI),
          angles[1] / 180 * Scalar(M_PI),
          angles[2] / 180 * Scalar(M_PI));
        Matrix33 R = euler_angles.template ToOrthoRotMat<2, 1, -3, 1>();
        R.transposeInPlace();
        R.col(2) *= -1;
        R.col(1) *= -1;
        extrinsic_params_set[id].rotation() = R;
      }
    }
  }

  Scalar focal_length_in_metre() const
  {
    return focal_length_in_metre_;
  }

  size_t number_of_strips() const
  {
    return number_of_strips_;
  }

  size_t number_of_cameras_in_strip() const
  {
    return number_of_cameras_in_strip_;
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

  Scalar pixel_size() const
  {
    return pixel_size_;
  }

  size_t number_of_points() const
  {
    return number_of_points_;
  }

  Scalar lateral_overlap_ratio() const
  {
    return lateral_overlap_ratio_;
  }

  Scalar longitudinal_overlap_ratio() const
  {
    return longitudinal_overlap_ratio_;
  }

  Scalar scene_max_height() const
  {
    return scene_max_height_;
  }

  Scalar camera_height_stddev() const
  {
    return camera_height_stddev_;
  }

  Scalar camera_planar_stddev() const
  {
    return camera_planar_stddev_;
  }

  Scalar camera_rotation_stddev() const
  {
    return camera_rotation_stddev_;
  }

private:
  /**
   *  相机的焦距米为单位
   */
  Scalar focal_length_in_metre_;
  /**
   *  航带数
   */
  size_t number_of_strips_;
  /**
   *  每条航带包含的相机数
   */
  size_t number_of_cameras_in_strip_;
  /**
   *  影像的地面分辨率
   */
  Scalar ground_resolution_;
  /**
   *  影像的宽度，像素为单位
   */
  ImageDimension image_width_;
  /**
   *  影像的高度，像素为单位
   */
  ImageDimension image_height_;
  /**
   *  像素大小，米为单位
   */
  Scalar pixel_size_;
  /**
   *  场景中空间点的数量
   */
  size_t number_of_points_;
  /**
   *  旁向重叠率
   */
  Scalar lateral_overlap_ratio_;
  /**
   *  航向重叠率
   */
  Scalar longitudinal_overlap_ratio_;
  /**
   *  场景中最高点的高度
   */
  Scalar scene_max_height_;
  /**
   *  相机高度分布标准差
   */
  Scalar camera_height_stddev_;
  /**
   *  相机平面方向分布标准差
   */
  Scalar camera_planar_stddev_;
  /**
   *  相机旋转角标准差，度为单位
   */
  Scalar camera_rotation_stddev_;
};

}
}
}

#endif
