#ifndef _HS_SFM_UTILITY_SYNTHETIC_SCENE_GENERATOR_HPP_
#define _HS_SFM_UTILITY_SYNTHETIC_SCENE_GENERATOR_HPP_

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_math/geometry/euler_angles.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/sfm_file_io.hpp"

namespace hs
{
namespace sfm
{

template <typename _Scalar, typename _ImageDimension>
class SceneGenerator
{
public:
typedef _Scalar Scalar;
typedef _ImageDimension ImageDimension;
typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
typedef typename ExtrinsicParams::Position Position;
typedef typename ExtrinsicParams::Rotation Rotation;
typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
typedef ImageParams<ImageDimension> Image;
typedef EIGEN_STD_VECTOR(Image) ImageContainer;
typedef EIGEN_VECTOR(Scalar, 3) Point3D;
typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;
typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

typedef int Err;

SceneGenerator(Scalar focal_length_in_metre,
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
               Scalar camera_rotation_stddev,
               Scalar north_west_angle)
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
       camera_rotation_stddev_(camera_rotation_stddev),
       north_west_angle_(north_west_angle) {}

Scalar GetFocalLengthInPixel() const
{
  return focal_length_in_metre_ / pixel_size_;
}

Err operator ()(IntrinsicParamsContainer& intrinsic_params_set,
                ExtrinsicParamsContainer& extrinsic_params_set,
                ImageContainer& images,
                Point3DContainer& points) const
{
  GenerateSceneCameras(intrinsic_params_set, extrinsic_params_set, images);
  GenerateScenePoints(number_of_points_, points);

  return 0;
}

inline Scalar GetFlightHeight() const
{
  return focal_length_in_metre_ * ground_resolution_ / pixel_size_;
}

void GetSceneDimensions(Scalar& scene_x_dimension,
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

void GenerateScenePoints(size_t number_of_points,
                         Point3DContainer& points) const
{
  points.resize(number_of_points);
  Scalar scene_x_dimension, scene_y_dimension, scene_z_dimension;
  GetSceneDimensions(scene_x_dimension, scene_y_dimension, scene_z_dimension);

  Scalar north_west_angle_radian = north_west_angle_ / 180 * Scalar(M_PI);
  Matrix33 north_west_rotaion;
  north_west_rotaion
    << cos(north_west_angle_radian), -sin(north_west_angle_radian), 0,
       sin(north_west_angle_radian),  cos(north_west_angle_radian), 0,
       0, 0, 1;

  Point3D max_point;
  max_point << scene_x_dimension * 0.5,
               scene_y_dimension * 0.5,
               scene_z_dimension;
  Point3D min_point;
  min_point << -scene_z_dimension * 0.5,
               -scene_y_dimension * 0.5,
               0;

  for (size_t i = 0; i < number_of_points; i++)
  {
    hs::math::random::UniformRandomVar<Scalar, 3>::Generate(
      min_point, max_point, points[i]);
    points[i] = north_west_rotaion * points[i];
  }
}

void GenerateScenePlanarPoints(size_t number_of_points,
                               Scalar height,
                               Point3DContainer& points) const
{
  points.resize(number_of_points);
  Scalar scene_x_dimension, scene_y_dimension, scene_z_dimension;
  GetSceneDimensions(scene_x_dimension, scene_y_dimension, scene_z_dimension);

  Scalar north_west_angle_radian = north_west_angle_ / 180 * Scalar(M_PI);
  Matrix33 north_west_rotaion;
  north_west_rotaion
    << cos(north_west_angle_radian), -sin(north_west_angle_radian), 0,
       sin(north_west_angle_radian),  cos(north_west_angle_radian), 0,
       0, 0, 1;

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
    points[i] = north_west_rotaion * points[i];
  }
}

void GenerateSceneCameras(IntrinsicParamsContainer& intrinsic_params_set,
                          ExtrinsicParamsContainer& extrinsic_params_set,
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

  intrinsic_params_set.resize(number_of_cameras,
                              IntrinsicParams(focal_length_in_metre_ /
                                              pixel_size_));
  extrinsic_params_set.resize(number_of_cameras);
  images.resize(number_of_cameras);

  Scalar scene_x_dimension, scene_y_dimension, scene_z_dimension;
  GetSceneDimensions(scene_x_dimension, scene_y_dimension, scene_z_dimension);

  Matrix33 camera_position_covariance = Matrix33::Identity();
  camera_position_covariance(0, 0) = camera_planar_stddev_;
  camera_position_covariance(1, 1) = camera_planar_stddev_;
  camera_position_covariance(2, 2) = camera_height_stddev_;
  Matrix33 camera_rotation_covariance = Matrix33::Identity();
  camera_rotation_covariance *= camera_rotation_stddev_;
  Scalar north_west_angle_radian = north_west_angle_ / 180 * Scalar(M_PI);
  Matrix33 north_west_rotation;
  north_west_rotation
    << cos(north_west_angle_radian), -sin(north_west_angle_radian), 0,
       sin(north_west_angle_radian),  cos(north_west_angle_radian), 0,
       0, 0, 1;
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
      extrinsic_params_set[id].position() =
        north_west_rotation * extrinsic_params_set[id].position();
      //extrins[id].m_c /= scale;
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
      extrinsic_params_set[id].rotation() =
        R * north_west_rotation.transpose();
    }
  }
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
/**
 *  飞行方向的北偏东角，度为单位
 */
Scalar north_west_angle_;
};

template <typename _Scalar, typename _ImageDimension>
class KeysGenerator
{
public:
typedef _Scalar Scalar;
typedef _ImageDimension ImageDimension;

typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinContainer;
typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
typedef typename ExtrinsicParams::Position Position;
typedef typename ExtrinsicParams::Rotation Rotation;
typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinContainer;
typedef EIGEN_VECTOR(Scalar, 3) Point3D;
typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;
typedef EIGEN_VECTOR(Scalar, 2) Point2D;
typedef EIGEN_STD_VECTOR(Point2D) Point2DSet;
typedef EIGEN_STD_VECTOR(Point2DSet) Point2DSetContainer;
typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
typedef CameraFunctions<Scalar> Camera;
typedef typename Camera::ProjectionMatrix ProjectionMatrix;
typedef ImageKeys<Scalar> Keys;
typedef EIGEN_STD_VECTOR(Keys) KeysContainer;

typedef int Err;

KeysGenerator(ImageDimension image_width, ImageDimension image_height)
  : image_width_(image_width), image_height_(image_height) {}

Err operator()(const IntrinContainer& intrinsic_params_set, 
               const ExtrinContainer& extrinsic_params_set, 
               const Point3DContainer& point,
               KeysContainer& keys_set, 
               TrackContainer& tracks,
               CameraViewContainer& camera_views) const
{
  size_t number_of_cameras = intrinsic_params_set.size();
  if (number_of_cameras != extrinsic_params_set.size())
  {
  return -1;
  }

  Point2DSetContainer key_sets(number_of_cameras);
  size_t number_of_points = point.size();
  tracks.resize(number_of_points);
  camera_views.resize(number_of_cameras);
  for (size_t i = 0; i < number_of_points; i++)
  {
    const Point3D& pt = point[i];
    for (size_t j = 0; j < number_of_cameras; j++)
    {
      ProjectionMatrix P =
        Camera::GetProjectionMatrix(intrinsic_params_set[j],
                                    extrinsic_params_set[j]);
      Point3D key_homogeneous = P.block(0, 0, 3, 3) * point[i] + 
                                P.block(0, 3, 3, 1);
      key_homogeneous /= key_homogeneous(2);
      if (key_homogeneous(0) > (-Scalar(image_width_) / 2) && 
          key_homogeneous(0) < ( Scalar(image_width_) / 2) &&
          key_homogeneous(1) > (-Scalar(image_height_) / 2) && 
          key_homogeneous(1) < ( Scalar(image_height_) / 2))
      {
        tracks[i].push_back(
          std::make_pair(j, key_sets[j].size()));
        camera_views[j].push_back(
          std::make_pair(i, key_sets[j].size()));
        key_sets[j].push_back(key_homogeneous.segment(0, 2));
      }
    }
  }
  keys_set.clear();
  for (size_t i = 0; i < number_of_cameras; i++)
  {
    keys_set.push_back(Keys(key_sets[i]));
  }

  return 0;
}

private:
ImageDimension image_width_;
ImageDimension image_height_;
};

class MatchGenerator
{
public:
typedef std::vector<std::pair<size_t, size_t> > Track;
typedef std::vector<Track> TrackContainer;
typedef int Err;

Err operator()(const TrackContainer& tracks,
               hs::sfm::MatchContainer& matches) const
{
  using namespace hs::sfm;
  size_t number_of_tracks = tracks.size();
  matches.clear();
  for (size_t i = 0; i < number_of_tracks; i++)
  {
    size_t number_of_views = tracks[i].size();
    for (size_t j = 0; j < number_of_views; j++)
    {
      for (size_t k = j + 1; k < number_of_views; k++)
      {
      ImagePair image_pair(tracks[i][j].first, tracks[i][k].first);
    
      matches[image_pair].push_back(
        KeyPair(tracks[i][j].second, tracks[i][k].second));
      }
    }
  }

  return 0;
}
};

template <typename _Scalar>
struct NoiseImageKeys
{
typedef _Scalar Scalar;
typedef ImageKeys<Scalar> ImgKeys;
typedef EIGEN_VECTOR(Scalar, 2) Key;
typedef EIGEN_MATRIX(Scalar, 2, 2) Covariance;
typedef int Err;

NoiseImageKeys(const Covariance& covariance)
  : covariance_(covariance) {}

Err operator()(ImgKeys& image_keys) const
{
  size_t number_of_image_keys = image_keys.size();
  for (size_t i = 0; i < number_of_image_keys; i++)
  {
  Key mean = image_keys[i];
  hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
    mean, covariance_, image_keys[i]);
  }

  return 0;
}

Covariance covariance_;
};

template <typename _Scalar>
struct NoiseExtrin
{
typedef _Scalar Scalar;
typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
typedef typename ExtrinsicParams::Position Position;
typedef typename ExtrinsicParams::Rotation Rotation;
typedef hs::math::geometry::EulerAngles<Scalar> EulerAnglesRotaion;
typedef EIGEN_MATRIX(Scalar, 3, 3) RotationMatrix;
typedef EIGEN_MATRIX(Scalar, 3, 3) Covariance;
typedef int Err;

NoiseExtrin(const Covariance& rotation_covariance,
            const Covariance& position_covariance)
  : rotation_covariance_(rotation_covariance),
    position_covariance_(position_covariance) {}

Err operator()(ExtrinsicParams& extrinsic_params) const
{
  Position mean_position = extrinsic_params.position();
  hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
  mean_position, position_covariance_, extrinsic_params.position());

  RotationMatrix R = extrinsic_params.rotation();
  R.transposeInPlace();
  EulerAnglesRotaion euler_angles;
  euler_angles.template FromOrthoRotMat<2, 1, -3, 1>(R);
  Position mean_rotation;
  mean_rotation << euler_angles[0],
  euler_angles[1],
  euler_angles[2];
  Position noised_angles;
  hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
  mean_rotation, rotation_covariance_, noised_angles);
  euler_angles[0] = noised_angles[0];
  euler_angles[1] = noised_angles[1];
  euler_angles[2] = noised_angles[2];
  R = euler_angles.template ToOrthoRotMat<2, 1, -3, 1>();
  R.transposeInPlace();
  extrinsic_params.rotation() = R;

  return 0;
}

Covariance rotation_covariance_;
Covariance position_covariance_;
};

}//sfm
}//hs

#endif
