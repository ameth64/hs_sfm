#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_VECTOR_FUNCTION_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_VECTOR_FUNCTION_HPP_

#include <vector>
#include <utility>

#include "hs_math/linear_algebra/eigen_macro.hpp"
//#include "hs_math/geometry/rotation.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

/**
 *  最简单的BA函数。将相机外参数与三维点作为X，特征点作为Y。
 *  三维点与特征点均为非齐次坐标，不考虑相机内参数。
 *  相机旋转使用轴角旋转。
 */
template <typename _Scalar>
class BANaiveVectorFunction
{
public:
  typedef _Scalar Scalar;

  /**
   *  X向量包含相机参数与三维点坐标，相机参数排在前边
   */
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) XVector;
  /**
   *  Y向量包含二维特征点坐标
   */
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) YVector;

  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;

  typedef typename XVector::Index Index;

  typedef int Err;

  /**
   *  每台相机的参数个数为6，其中旋转占3个，位置占3个
   */
  static const Index params_per_camera_ = 6;

  /**
   *  三维点参数个数为3
   */
  static const Index params_per_point_ = 3;

  /**
   *  特征点参数个数为2
   */
  static const Index params_per_feature_ = 2;

  /**
   *  表示特征点与相机、三维点之间的映射关系
   *  pair的first值表示相机序号，second值表示三维点序号
   */
  typedef std::pair<Index, Index> FeatureMap;
  typedef std::vector<FeatureMap> FeatureMapContainer;

public:
  BANaiveVectorFunction() {}
  BANaiveVectorFunction(Index number_of_cameras,
                        Index number_of_points_,
                        Index number_of_features,
                        const FeatureMapContainer& feature_maps)
    : number_of_cameras_(number_of_cameras),
      number_of_points_(number_of_points_),
      number_of_features_(number_of_features),
      feature_maps_(feature_maps) {}

  inline Index number_of_cameras() const
  {
    return number_of_cameras_;
  }

  inline void set_number_of_cameras(Index number_of_cameras)
  {
    number_of_cameras_ = number_of_cameras;
  }

  inline Index number_of_points() const
  {
    return number_of_points_;
  }

  inline void set_number_of_points(Index number_of_points)
  {
    number_of_points_ = number_of_points;
  }

  inline Index number_of_features() const
  {
    return number_of_features_;
  }

  inline void set_number_of_features(Index number_of_features)
  {
    number_of_features_ = number_of_features;
  }

  inline const FeatureMapContainer& feature_maps() const
  {
    return feature_maps_;
  }

  inline void set_feature_maps(const FeatureMapContainer& feature_maps)
  {
    feature_maps_ = feature_maps;
  }

  inline Index GetXSize() const
  {
    return GetCameraParamsSize() + GetPointParamsSize();
  }

  inline Index GetYSize() const
  {
    return number_of_features_ * params_per_feature_;
  }

  inline Index GetCameraParamsSize() const
  {
    return number_of_cameras_ * params_per_camera_;
  }

  inline Index GetPointParamsSize() const
  {
    return number_of_points_ * params_per_point_;
  }

  Err operator() (const XVector& x, YVector& y) const
  {
    Index x_size = x.rows();
    if (x_size != GetXSize()) return -1;
    if (number_of_features_ != Index(feature_maps_.size())) return -1;
    Index y_size = GetYSize();
    y.resize(y_size);

    Index camera_params_size = GetCameraParamsSize();

    for (Index i = 0; i < number_of_features_; i++)
    {
      Index camera_id = feature_maps_[i].first;
      Index point_id = feature_maps_[i].second;

      Vector3 p = x.segment(camera_params_size + point_id * params_per_point_,
                            params_per_point_);

      Vector3 r = x.segment(camera_id * params_per_camera_, 3);
      Vector3 t = x.segment(camera_id * params_per_camera_ + 3, 3);

      Vector2 feat;
      PointProjectToFeature(r, t, p, feat);

      y.segment(i * params_per_feature_, params_per_feature_) = feat;
    }

    return 0;
  }

  inline static Err PointProjectToFeature(const Vector3& rotation,
                                          const Vector3& translation,
                                          const Vector3& point,
                                          Vector2& feature)
  {
    //hs::math::geometry::Rot3D<Scalar> R(r);
    //Vec3 ptCam = R * p + t;

    //轴角旋转加上平移
    Scalar theta = rotation.norm();
    if (theta == Scalar(0))
    {
      Vector3 feature_homogeneous = point + translation;
      feature_homogeneous /= feature_homogeneous[2];
      feature = feature_homogeneous.segment(0, 2);
    }
    else
    {
      Vector3 normalized_rotation = rotation / theta;
      Vector3 feature_homogeneous =
        cos(theta) * point +
        sin(theta) * normalized_rotation.cross(point) +
        (1 - cos(theta)) * point.dot(normalized_rotation) *
        normalized_rotation + translation;
      feature_homogeneous /= feature_homogeneous[2];
      feature = feature_homogeneous.segment(0, 2);
    }

    return 0;
  }

private:
  Index number_of_cameras_;
  Index number_of_points_;
  Index number_of_features_;
  FeatureMapContainer feature_maps_;
};

}//ba
}//sfm
}//hs

#endif
