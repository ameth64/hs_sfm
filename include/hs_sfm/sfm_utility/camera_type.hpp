#ifndef _HS_SFM_UTILITY_CAMERA_TYPE_HPP_
#define _HS_SFM_UTILITY_CAMERA_TYPE_HPP_

#include <bitset>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_math/geometry/rotation.hpp"

namespace hs
{
namespace sfm
{

template <typename _Scalar>
class CameraIntrinsicParams
{
public:
  typedef _Scalar Scalar;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  /**
   *  内参数矩阵，表示相机内部的固有参数
   */
  typedef Matrix33 KMatrix;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  CameraIntrinsicParams()
    : focal_length_(0),
      skew_(0),
      principal_point_x_(0),
      principal_point_y_(0),
      pixel_ratio_(1) {}
  CameraIntrinsicParams(Scalar focal_length,
        Scalar skew = 0,
        Scalar principal_point_x = 0,
        Scalar principal_point_y = 0,
        Scalar pixel_ratio = 1)
    : focal_length_(focal_length),
      skew_(skew),
      principal_point_x_(principal_point_x),
      principal_point_y_(principal_point_y),
      pixel_ratio_(pixel_ratio) {}
  CameraIntrinsicParams(const KMatrix& K)
    : focal_length_(K(0, 0)), 
      skew_(K(0, 1)),
      principal_point_x_(K(0, 2)),
      principal_point_y_(K(1, 2)),
      pixel_ratio_(K(0, 0) / K(1, 1)) {}

  inline KMatrix GetKMatrix() const
  {
    KMatrix K;
    K << focal_length_, skew_, principal_point_x_,
       0, focal_length_ / pixel_ratio_, principal_point_y_,
       0, 0, 1;

    return K;
  }

  inline operator KMatrix() const
  {
    return GetKMatrix();
  }

  Scalar focal_length() const
  {
    return focal_length_;
  }

  void set_focal_length(Scalar focal_length)
  {
    focal_length_ = focal_length;
  }

  Scalar skew() const
  {
    return skew_;
  }

  void set_skew(Scalar skew)
  {
    skew_ = skew;
  }

  Scalar principal_point_x() const
  {
    return principal_point_x_;
  }

  void set_principal_point_x(Scalar principal_point_x)
  {
    principal_point_x_ = principal_point_x;
  }

  Scalar principal_point_y() const
  {
    return principal_point_y_;
  }

  void set_principal_point_y(Scalar principal_point_y)
  {
    principal_point_y_ = principal_point_y;
  }

  Scalar pixel_ratio() const
  {
    return pixel_ratio_;
  }

  void set_pixel_ratio(Scalar pixel_ratio)
  {
    pixel_ratio_ = pixel_ratio;
  }

private:
  /**
   *  焦距
   */
  Scalar focal_length_;
  /**
   *  扭曲参数
   */
  Scalar skew_;
  /**
   *  主点偏移
   */
  Scalar principal_point_x_;
  Scalar principal_point_y_;
  /**
   *  x方向与y方向上像素量纲比值
   */
  Scalar pixel_ratio_;
};

template <typename _Scalar>
class CameraExtrinsicParams
{
public:
  typedef _Scalar Scalar;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_MATRIX(Scalar, 3, 4) Matrix34;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef hs::math::geometry::Rotation3D<Scalar> Rotation;
  typedef Vector3 Position;
  /**
   *  外参数矩阵,表示相机运动
   */
  typedef Matrix34 MotionMatrix;

  CameraExtrinsicParams() {rotation_.SetIdentity();position_.setZero();}
  CameraExtrinsicParams(const Rotation& rotation, const Position& position)
    : rotation_(rotation), position_(position) {}

  inline MotionMatrix GetMotionMatrix() const
  {
    MotionMatrix M;
    M.block(0, 0, 3, 3) = Matrix33(rotation_);
    M.block(0, 3, 3, 1) = -Matrix33(rotation_) * position_;
    return M;
  }

  inline operator MotionMatrix() const
  {
    return GetMotionMatrix();
  }

  const Rotation& rotation() const
  {
    return rotation_;
  }

  Rotation& rotation()
  {
    return rotation_;
  }

  const Position& position() const
  {
    return position_;
  }

  Position& position()
  {
    return position_;
  }

private:
  /**
   *  相机朝向。
   *  表示由世界坐标系到相机坐标系的旋转。
   */
  Rotation rotation_;
  /**
   *  相机位置。
   *  表示相机在世界坐标系下的位置。
   */
  Position position_;
};

template <typename _Dimension>
struct ImageParams
{
  typedef _Dimension Dimension;
  
  size_t m_id;
  Dimension m_width;
  Dimension m_height;

  std::string m_path;
};

template <typename _Scalar>
class CameraFunctions
{
public:
  typedef _Scalar Scalar;

  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;

  typedef EIGEN_MATRIX(Scalar, 3, 4) ProjectionMatrix;
  typedef EIGEN_MATRIX(Scalar, 3, 3) RotationMatrix;

  static ProjectionMatrix GetProjectionMatrix(
    const IntrinsicParams& intrinsic_params,
    const ExtrinsicParams& extrinsic_params)
  {
    ProjectionMatrix projection_matrix = extrinsic_params.GetMotionMatrix();
    return intrinsic_params.GetKMatrix() * projection_matrix;
  }
};

}
}

#endif
