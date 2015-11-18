#ifndef _HS_SFM_SFM_UTILITY_CAMERA_ROTATED_BOUNDING_BOX_INTERSECTOR_HPP_
#define _HS_SFM_SFM_UTILITY_CAMERA_ROTATED_BOUNDING_BOX_INTERSECTOR_HPP_

#include "hs_math/geometry/plane_from_3points.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/projective_functions.hpp"

namespace hs
{
namespace sfm
{

/**
 *  计算相机的投影四棱锥与经过旋转的立方体矩形框是否相交的模板类。
 *
 *  计算方法如下：
 *  - 计算立方体矩形框是否存在投影四棱锥内的顶点，若存在，则相交；
 *  - 计算相机的投影四棱锥的四条边与立方体矩形框的交点，
 *    若存在交点在立方体矩形框内，则相交。
 */
template <typename _Scalar>
class CameraRotatedBoundingBoxIntersector
{
public:
  typedef _Scalar Scalar;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;

private:
  typedef hs::sfm::ProjectiveFunctions<Scalar> ProjectiveFunctions;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef EIGEN_VECTOR(Scalar, 4) Vector4;
  typedef EIGEN_MATRIX(Scalar, 4, 4) Matrix44;
  typedef hs::math::geometry::PlaneFrom3Points<Scalar> PlaneFrom3Points;

public:
  bool operator()(const IntrinsicParams& intrinsic_params,
                  const ExtrinsicParams& extrinsic_params,
                  size_t width, size_t height,
                  const Point& center,
                  const Point& min, const Point& max,
                  const RMatrix& rmatrix) const
  {
    Point corners[8] =
    {
      Point(min[0], min[1], min[2]),
      Point(min[0], min[1], max[2]),
      Point(min[0], max[1], min[2]),
      Point(min[0], max[1], max[2]),
      Point(max[0], min[1], min[2]),
      Point(max[0], min[1], max[2]),
      Point(max[0], max[1], min[2]),
      Point(max[0], max[1], max[2])
    };

    bool is_all_corners_outside = true;
    for (auto& corner : corners)
    {
      corner = rmatrix.transpose() * corner + center;
      Key key;
      key = ProjectiveFunctions::WorldPointProjectToImageKey(
        intrinsic_params, extrinsic_params, corner);
      if (key[0] >= Scalar(0) && key[0] < Scalar(width) &&
          key[1] >= Scalar(0) && key[1] < Scalar(height))
      {
        key = ProjectiveFunctions::WorldPointProjectToImageKeyNoDistort(
          intrinsic_params, extrinsic_params, corner);
        if (key[0] >= Scalar(0) && key[0] < Scalar(width) &&
            key[1] >= Scalar(0) && key[1] < Scalar(height))
        {
          is_all_corners_outside = false;
          break;
        }
      }
    }
    if (is_all_corners_outside)
    {
      Vector4 position;
      position.template segment<3>(0) = extrinsic_params.position();
      position[3] = Scalar(1);
      Vector4 camera_corners[4] =
      {
        Vector4(-Scalar(width) / Scalar(2),
              -Scalar(height) / Scalar(2),
              Scalar(intrinsic_params.focal_length()),
              Scalar(1)),
        Vector4(Scalar(width) / Scalar(2),
              -Scalar(height) / Scalar(2),
              Scalar(intrinsic_params.focal_length()),
              Scalar(1)),
        Vector4(Scalar(width) / Scalar(2),
              Scalar(height) / Scalar(2),
              Scalar(intrinsic_params.focal_length()),
              Scalar(1)),
        Vector4(-Scalar(width) / Scalar(2),
              Scalar(height) / Scalar(2),
              Scalar(intrinsic_params.focal_length()),
              Scalar(1))
      };
      RMatrix rmatrix_camera = extrinsic_params.rotation();
      Point translate = -rmatrix_camera * extrinsic_params.position();
      for (int i = 0; i < 4; i++)
      {
        camera_corners[i].template segment<3>(0) =
          rmatrix_camera.transpose() *
          (camera_corners[i].template segment<3>(0) - translate);
      }
      Vector4 plane_left =
        PlaneFrom3Points()(corners[0], corners[1], corners[2]);
      Vector4 plane_right =
        PlaneFrom3Points()(corners[4], corners[5], corners[6]);
      Vector4 plane_bottom =
        PlaneFrom3Points()(corners[0], corners[1], corners[5]);
      Vector4 plane_top =
        PlaneFrom3Points()(corners[2], corners[3], corners[7]);
      Vector4 plane_back =
        PlaneFrom3Points()(corners[2], corners[4], corners[6]);
      Vector4 plane_front =
        PlaneFrom3Points()(corners[3], corners[5], corners[7]);

      for (int i = 0; i < 4; i++)
      {
        Matrix44 line = position * camera_corners[i].transpose() -
                        camera_corners[i] * position.transpose();
        for (int j = 0; j < 6; j++)
        {
          Vector4 point_left = line * plane_left;
          if (std::abs(point_left[3]) > 1e-10)
          {
            point_left /= point_left[3];
            Point point_left_r = point_left.template segment<3>(0);
            point_left_r = rmatrix * (point_left_r - center);
            if (point_left_r[1] > min[1] && point_left_r[1] < max[1] &&
                point_left_r[2] > min[2] && point_left_r[2] < max[2])
            {
              return true;
            }
          }
          Vector4 point_right = line * plane_right;
          if (std::abs(point_right[3]) > 1e-10)
          {
            point_right /= point_right[3];
            Point point_right_r = point_right.template segment<3>(0);
            point_right_r = rmatrix * (point_right_r - center);
            if (point_right_r[1] > min[1] && point_right_r[1] < max[1] &&
                point_right_r[2] > min[2] && point_right_r[2] < max[2])
            {
              return true;
            }
          }
          Vector4 point_bottom = line * plane_bottom;
          if (std::abs(point_bottom[3]) > 1e-10)
          {
            point_bottom /= point_bottom[3];
            Point point_bottom_r = point_bottom.template segment<3>(0);
            point_bottom_r = rmatrix * (point_bottom_r - center);
            if (point_bottom_r[0] > min[0] && point_bottom_r[0] < max[0] &&
                point_bottom_r[2] > min[2] && point_bottom_r[2] < max[2])
            {
              return true;
            }
          }
          Vector4 point_top = line * plane_top;
          if (std::abs(point_top[3]) > 1e-10)
          {
            point_top /= point_top[3];
            Point point_top_r = point_top.template segment<3>(0);
            point_top_r = rmatrix * (point_top_r - center);
            if (point_top_r[0] > min[0] && point_top_r[0] < max[0] &&
                point_top_r[2] > min[2] && point_top_r[2] < max[2])
            {
              return true;
            }
          }
          Vector4 point_back = line * plane_back;
          if (std::abs(point_back[3]) > 1e-10)
          {
            point_back /= point_back[3];
            Point point_back_r = point_back.template segment<3>(0);
            point_back_r = rmatrix * (point_back_r - center);
            if (point_back_r[0] > min[0] && point_back_r[0] < max[0] &&
                point_back_r[1] > min[1] && point_back_r[1] < max[1])
            {
              return true;
            }
          }
          Vector4 point_front = line * plane_front;
          if (std::abs(point_front[3]) > 1e-10)
          {
            point_front /= point_front[3];
            Point point_front_r = point_front.template segment<3>(0);
            point_front_r = rmatrix * (point_front_r - center);
            if (point_front_r[0] > min[0] && point_front_r[0] < max[0] &&
                point_front_r[1] > min[1] && point_front_r[1] < max[1])
            {
              return true;
            }
          }
        }
      }

      return false;
    }
    else
    {
      return true;
    }
  }
};

}
}

#endif
