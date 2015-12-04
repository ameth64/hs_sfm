#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_VECTOR_FUNCTION_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_VECTOR_FUNCTION_HPP_

#include <vector>
#include <utility>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/sfm_utility/radial_distortor.hpp"
#include "hs_sfm/sfm_utility/decentering_distortor.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_common_types.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedVectorFunction
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) XVector;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) YVector;

  typedef typename XVector::Index Index;

  static const Index extrinsic_params_per_image_ = 6;
  static const Index params_per_point_ = 3;
  static const Index params_per_key_ = 2;

  typedef std::pair<Index, Index> FeatureMap;
  typedef std::vector<FeatureMap> FeatureMapContainer;

  typedef std::vector<Index> ImageCameraMap;

  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) VectorX;

  typedef EIGEN_VECTOR(Scalar, params_per_point_) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;
  typedef EIGEN_VECTOR(Scalar, extrinsic_params_per_image_) Image;
  typedef EIGEN_STD_VECTOR(Image) ImageContainer;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) Camera;
  typedef EIGEN_STD_VECTOR(Camera) CameraContainer;

public:
  CameraSharedVectorFunction()
    : number_of_images_(0),
      number_of_points_(0),
      number_of_keys_(0),
      number_of_cameras_(0),
      feature_maps_(),
      image_camera_map_(),
      intrinsic_computations_mask_(),
      point_constraints_(),
      image_constraints_(),
      camera_constraints_() {}

  inline Index number_of_cameras() const
  {
    return number_of_cameras_;
  }
  inline void set_number_of_cameras(Index number_of_cameras)
  {
    number_of_cameras_ = number_of_cameras;
  }

  inline Index number_of_images() const
  {
    return number_of_images_;
  }

  inline void set_number_of_images(Index number_of_images)
  {
    number_of_images_ = number_of_images;
  }

  inline Index number_of_points() const
  {
    return number_of_points_;
  }
  inline void set_number_of_points(Index number_of_points)
  {
    number_of_points_ = number_of_points;
  }

  inline Index number_of_keys() const
  {
    return number_of_keys_;
  }
  inline void set_number_of_keys(Index number_of_keys)
  {
    number_of_keys_ = number_of_keys;
  }

  inline const FeatureMapContainer& feature_maps() const
  {
    return feature_maps_;
  }
  inline void set_feature_maps(const FeatureMapContainer& feature_maps)
  {
    feature_maps_ = feature_maps;
  }

  inline const ImageCameraMap& image_camera_map() const
  {
    return image_camera_map_;
  }
  inline void set_image_camera_map(const ImageCameraMap& image_camera_map)
  {
    image_camera_map_ = image_camera_map;
  }

  inline const IntrinsicComputationsMask& intrinsic_computations_mask() const
  {
    return intrinsic_computations_mask_;
  }
  inline IntrinsicComputationsMask& intrinsic_computations_mask()
  {
    return intrinsic_computations_mask_;
  }

  inline const PointConstraintContainer& point_constraints() const
  {
    return point_constraints_;
  }
  inline PointConstraintContainer& point_constraints()
  {
    return point_constraints_;
  }

  inline const ImageConstraintContainer& image_constraints() const
  {
    return image_constraints_;
  }
  inline ImageConstraintContainer& image_constraints()
  {
    return image_constraints_;
  }

  inline const CameraConstraintContainer& camera_constraints() const
  {
    return camera_constraints_;
  }
  inline CameraConstraintContainer& camera_constraints()
  {
    return camera_constraints_;
  }

  inline Err set_fix_points(const PointContainer& fix_points)
  {
    if (Index(fix_points.size()) != number_of_points_)
    {
      return -1;
    }
    fix_mask_.set(FIX_POINTS);
    fix_points_ = fix_points;
    return 0;
  }

  inline void unset_fix_points()
  {
    fix_mask_.reset(FIX_POINTS);
    PointContainer().swap(fix_points_);
  }

  inline bool is_fix_points() const
  {
    return fix_mask_[FIX_POINTS];
  }

  inline const PointContainer& fix_points() const
  {
    return fix_points_;
  }

  inline Err set_fix_images(const ImageContainer& fix_images)
  {
    if (Index(fix_images.size()) != number_of_images_)
    {
      return -1;
    }
    fix_mask_.set(FIX_IMAGES);
    fix_images_ = fix_images;
    return 0;
  }

  inline void unset_fix_images()
  {
    fix_mask_.reset(FIX_IMAGES);
    ImageContainer().swap(fix_images_);
  }

  inline bool is_fix_images() const
  {
    return fix_mask_[FIX_IMAGES];
  }

  inline const ImageContainer& fix_images() const
  {
    return fix_images_;
  }

  inline Err set_fix_cameras(const CameraContainer& fix_cameras)
  {
    if (Index(fix_cameras.size()) != number_of_cameras_)
    {
      return -1;
    }
    fix_mask_.set(FIX_CAMERAS);
    fix_cameras_ = fix_cameras;
    return 0;
  }

  inline void unset_fix_cameras()
  {
    fix_mask_.reset(FIX_CAMERAS);
    CameraContainer().swap(fix_cameras_);
  }

  inline bool is_fix_cameras() const
  {
    return fix_mask_[FIX_CAMERAS];
  }

  inline const CameraContainer& fix_cameras() const
  {
    return fix_cameras_;
  }

  inline Index GetExtrinsicParamsSize() const
  {
    return number_of_images_ * extrinsic_params_per_image_;
  }

  inline Index GetPointParamsSize() const
  {
    return number_of_points_ * params_per_point_;
  }

  inline Index GetIntrinsicParamsSizePerCamera() const
  {
    Index params_per_camera = 0;
    if (intrinsic_computations_mask_[COMPUTE_RADIAL_DISTORTION])
    {
      params_per_camera += 3;
    }
    if (intrinsic_computations_mask_[COMPUTE_DECENTERING_DISTORTION])
    {
      params_per_camera += 2;
    }
    if (intrinsic_computations_mask_[COMPUTE_INTRINSIC_PARAMS])
    {
      params_per_camera += 5;
    }
    return params_per_camera;
  }

  inline Index GetIntrinsicParamsSize() const
  {
    return number_of_cameras_ * GetIntrinsicParamsSizePerCamera();
  }

  inline Index GetYKeysSize() const
  {
    return number_of_keys_ * params_per_key_;
  }

  inline Index GetYPointConstraintsSize() const
  {
    if (!is_fix_points())
    {
      Index point_constraints_size = 0;
      auto itr_point = point_constraints_.begin();
      auto itr_point_end = point_constraints_.end();
      for (; itr_point != itr_point_end; ++itr_point)
      {
        point_constraints_size += itr_point->mask.count();
      }
      return point_constraints_size;
    }
    else
    {
      return 0;
    }
  }

  inline Index GetYImageConstraintsSize() const
  {
    if (!is_fix_images())
    {
      Index image_constraints_size = 0;
      auto itr_image = image_constraints_.begin();
      auto itr_image_end = image_constraints_.end();
      for (; itr_image != itr_image_end; ++itr_image)
      {
        if (itr_image->mask[IMAGE_CONSTRAIN_ROTATION])
        {
          image_constraints_size += 3;
        }
        if (itr_image->mask[IMAGE_CONSTRAIN_POSITION_X])
        {
          image_constraints_size++;
        }
        if (itr_image->mask[IMAGE_CONSTRAIN_POSITION_Y])
        {
          image_constraints_size++;
        }
        if (itr_image->mask[IMAGE_CONSTRAIN_POSITION_Z])
        {
          image_constraints_size++;
        }
      }
      return image_constraints_size;
    }
    else
    {
      return 0;
    }
  }

  inline Index GetYCameraConstraintsSize() const
  {
    if (!is_fix_cameras())
    {
      Index camera_constraints_size = 0;
      auto itr_camera = camera_constraints_.begin();
      auto itr_camera_end = camera_constraints_.end();
      for (; itr_camera != itr_camera_end; ++itr_camera)
      {
        camera_constraints_size += itr_camera->radial_mask.count() +
                                   itr_camera->decentering_mask.count() +
                                   itr_camera->intrinsic_mask.count();
      }
      return camera_constraints_size;
    }
    else
    {
      return 0;
    }
  }

  inline Index GetXSize() const
  {
    return ((is_fix_points() ? 0 : GetPointParamsSize()) +
            (is_fix_images() ? 0 : GetExtrinsicParamsSize()) +
            (is_fix_cameras() ? 0 : GetIntrinsicParamsSize()));
  }

  inline Index GetYSize() const
  {
    return (GetYKeysSize() +
            GetYPointConstraintsSize() +
            GetYImageConstraintsSize() +
            GetYCameraConstraintsSize());
  }

  Point GetPoint(Index point_id, const XVector& x) const
  {
    Point point;
    if (is_fix_points())
    {
      point = fix_points_[point_id];
    }
    else
    {
      point =
        x.template segment<params_per_point_>(point_id * params_per_point_);
    }
    return point;
  }

  Image GetImage(Index image_id, const XVector& x) const
  {
    Image image;
    if (is_fix_images())
    {
      image = fix_images_[image_id];
    }
    else
    {
      Index image_begin = is_fix_points() ? 0 : GetPointParamsSize();
      image = x.template segment<extrinsic_params_per_image_>(
                image_begin + image_id * extrinsic_params_per_image_);
    }
    return image;
  }

  Camera GetCamera(Index camera_id, const XVector& x) const
  {
    Index camera_size = GetIntrinsicParamsSizePerCamera();
    Camera camera(camera_size);
    if (is_fix_cameras())
    {
      camera = fix_cameras_[camera_id];
    }
    else
    {
      Index camera_begin = (is_fix_points() ? 0 : GetPointParamsSize()) +
                           (is_fix_images() ? 0 : GetExtrinsicParamsSize());
      camera = x.segment(camera_begin + camera_id * camera_size, camera_size);
    }
    return camera;
  }

  Err operator() (const XVector& x, YVector& y) const
  {
    Index x_size = x.rows();
    if (x_size != GetXSize()) return -1;
    if (number_of_keys_ != Index(feature_maps_.size())) return -1;
    if (number_of_images_ != Index(image_camera_map_.size())) return -1;

    Index y_size = GetYSize();
    y.resize(y_size);

    //keys segment
    Index point_params_size = GetPointParamsSize();
    Index extrinsic_params_size = GetExtrinsicParamsSize();
    Index intrinsic_params_size_per_camera =
      GetIntrinsicParamsSizePerCamera();
    Index x_point_begin = 0;
    Index x_image_begin = is_fix_points() ? 0 : point_params_size;
    Index x_camera_begin = (is_fix_points() ? 0 : point_params_size) +
                           (is_fix_images() ? 0 : extrinsic_params_size);
    for (Index i = 0; i < number_of_keys_; i++)
    {
      Index image_id = feature_maps_[i].first;
      Index point_id = feature_maps_[i].second;
      Index camera_id = image_camera_map_[image_id];

      Vector3 point;
      if (is_fix_points())
      {
        point = fix_points_[point_id];
      }
      else
      {
        point = x.segment(point_id * params_per_point_, params_per_point_);
      }
      Vector3 rotation, translation;
      if (is_fix_images())
      {
        rotation = fix_images_[image_id].segment(0, 3);
        translation = fix_images_[image_id].segment(3, 3);
      }
      else
      {
        rotation =
          x.segment(x_image_begin +
                    image_id * extrinsic_params_per_image_, 3);
        translation =
          x.segment(x_image_begin +
                    image_id * extrinsic_params_per_image_ + 3, 3);
      }

      VectorX intrinsic_params(intrinsic_params_size_per_camera);
      if (is_fix_cameras())
      {
        intrinsic_params = fix_cameras_[camera_id];
      }
      else
      {
        intrinsic_params =
          x.segment(x_camera_begin +
                    camera_id * intrinsic_params_size_per_camera,
                    intrinsic_params_size_per_camera);
      }

      Vector2 image_key = WorldPointToImageKey(point,
                                               rotation,
                                               translation,
                                               intrinsic_params);

      y.segment(i * params_per_key_, params_per_key_) = image_key;
    }

    Index y_offset = GetYKeysSize();

    if (!is_fix_points())
    {
      auto itr_point_constraint = point_constraints_.begin();
      auto itr_point_constraint_end = point_constraints_.end();
      for (; itr_point_constraint != itr_point_constraint_end;
           ++itr_point_constraint)
      {
        Index point_id = Index(itr_point_constraint->point_id);
        Index x_offset = x_point_begin + point_id * params_per_point_;
        if (itr_point_constraint->mask[POINT_CONSTRAIN_X])
        {
          y[y_offset] = x[x_offset + 0];
          y_offset++;
        }
        if (itr_point_constraint->mask[POINT_CONSTRAIN_Y])
        {
          y[y_offset] = x[x_offset + 1];
          y_offset++;
        }
        if (itr_point_constraint->mask[POINT_CONSTRAIN_Z])
        {
          y[y_offset] = x[x_offset + 2];
          y_offset++;
        }
      }
    }
    if (!is_fix_images())
    {
      auto itr_image_constraint = image_constraints_.begin();
      auto itr_image_constraint_end = image_constraints_.end();
      for (; itr_image_constraint != itr_image_constraint_end;
           ++itr_image_constraint)
      {
        Index image_id = Index(itr_image_constraint->image_id);
        Index x_offset = x_image_begin + image_id * extrinsic_params_per_image_;
        if (itr_image_constraint->mask[IMAGE_CONSTRAIN_ROTATION])
        {
          y.segment(y_offset, 3) = x.segment(x_offset, 3);
          y_offset += 3;
        }
        if (itr_image_constraint->mask[IMAGE_CONSTRAIN_POSITION_X])
        {
          y[y_offset] = x[x_offset + 3];
          y_offset++;
        }
        if (itr_image_constraint->mask[IMAGE_CONSTRAIN_POSITION_Y])
        {
          y[y_offset] = x[x_offset + 4];
          y_offset++;
        }
        if (itr_image_constraint->mask[IMAGE_CONSTRAIN_POSITION_Z])
        {
          y[y_offset] = x[x_offset + 5];
          y_offset++;
        }
      }
    }
    if (!is_fix_cameras())
    {
      auto itr_camera_constraint = camera_constraints_.begin();
      auto itr_camera_constraint_end = camera_constraints_.end();
      for (; itr_camera_constraint != itr_camera_constraint_end;
           ++itr_camera_constraint)
      {
        if ((itr_camera_constraint->radial_mask.any() &&
             !intrinsic_computations_mask_[COMPUTE_RADIAL_DISTORTION]) ||
            (itr_camera_constraint->decentering_mask.any() &&
             !intrinsic_computations_mask_[COMPUTE_DECENTERING_DISTORTION]) ||
            (itr_camera_constraint->intrinsic_mask.any() &&
             !intrinsic_computations_mask_[COMPUTE_INTRINSIC_PARAMS]))
        {
          return -1;
        }
        Index camera_id = Index(itr_camera_constraint->camera_id);
        Index x_offset = x_camera_begin +
                         camera_id * intrinsic_params_size_per_camera;
        if (itr_camera_constraint->radial_mask[RADIAL_CONSTRAIN_K1])
        {
          y[y_offset] = x[x_offset + 0];
          y_offset++;
        }
        if (itr_camera_constraint->radial_mask[RADIAL_CONSTRAIN_K2])
        {
          y[y_offset] = x[x_offset + 1];
          y_offset++;
        }
        if (itr_camera_constraint->radial_mask[RADIAL_CONSTRAIN_K3])
        {
          y[y_offset] = x[x_offset + 2];
          y_offset++;
        }

        if (intrinsic_computations_mask_[COMPUTE_RADIAL_DISTORTION])
        {
          x_offset += 3;
        }

        if (itr_camera_constraint->decentering_mask[DECENTERING_CONSTRAIN_D1])
        {
          y[y_offset] = x[x_offset + 0];
          y_offset++;
        }
        if (itr_camera_constraint->decentering_mask[DECENTERING_CONSTRAIN_D2])
        {
          y[y_offset] = x[x_offset + 1];
          y_offset++;
        }

        if (intrinsic_computations_mask_[COMPUTE_DECENTERING_DISTORTION])
        {
          x_offset += 2;
        }

        if (itr_camera_constraint->intrinsic_mask[
              INTRINSIC_CONSTRAIN_FOCAL_LENGTH])
        {
          y[y_offset] = x[x_offset + 0];
          y_offset++;
        }
        if (itr_camera_constraint->intrinsic_mask[
              INTRINSIC_CONSTRAIN_SKEW])
        {
          y[y_offset] = x[x_offset + 1];
          y_offset++;
        }
        if (itr_camera_constraint->intrinsic_mask[
              INTRINSIC_CONSTRAIN_PRINCIPAL_X])
        {
          y[y_offset] = x[x_offset + 2];
          y_offset++;
        }
        if (itr_camera_constraint->intrinsic_mask[
              INTRINSIC_CONSTRAIN_PRINCIPAL_Y])
        {
          y[y_offset] = x[x_offset + 3];
          y_offset++;
        }
        if (itr_camera_constraint->intrinsic_mask[
              INTRINSIC_CONSTRAIN_PIXEL_RATIO])
        {
          y[y_offset] = x[x_offset + 4];
          y_offset++;
        }
      }
    }

    return 0;
  }

  Vector2 WorldPointToImageKey(const Vector3& world_point,
                               const Vector3& rotation,
                               const Vector3& translation,
                               const VectorX& intrinsic_params) const
  {
    Vector2 normalized_key;
    WorldPointToNormalizedKey(world_point,
                              rotation,
                              translation,
                              normalized_key);

    Index offset = 0;
    Vector2 radial_delta = Vector2::Zero();
    if (intrinsic_computations_mask_[COMPUTE_RADIAL_DISTORTION])
    {
      Scalar k1 = intrinsic_params[offset + 0];
      Scalar k2 = intrinsic_params[offset + 1];
      Scalar k3 = intrinsic_params[offset + 2];
      GetRadialDistortDelta(k1, k2, k3, normalized_key,  radial_delta);
      offset += 3;
    }
    Vector2 decentering_delta = Vector2::Zero();
    if (intrinsic_computations_mask_[COMPUTE_DECENTERING_DISTORTION])
    {
      Scalar d1 = intrinsic_params[offset + 0];
      Scalar d2 = intrinsic_params[offset + 1];
      GetDecenteringDistortDelta(d1, d2, normalized_key, decentering_delta);
      offset += 2;
    }
    normalized_key += radial_delta + decentering_delta;
    Vector2 image_key;
    if (intrinsic_computations_mask_[COMPUTE_INTRINSIC_PARAMS])
    {
      Scalar focal_length = intrinsic_params[offset + 0];
      Scalar skew = intrinsic_params[offset + 1];
      Scalar principal_x = intrinsic_params[offset + 2];
      Scalar principal_y = intrinsic_params[offset + 3];
      Scalar pixel_ratio = intrinsic_params[offset + 4];
      NormalizedKeyToImageKey(focal_length,
                              skew,
                              principal_x,
                              principal_y,
                              pixel_ratio,
                              normalized_key,
                              image_key);
    }
    else
    {
      image_key = normalized_key;
    }

    return image_key;
  }

  inline static void WorldPointToNormalizedKey(const Vector3& world_point,
                                               const Vector3& rotation,
                                               const Vector3& translation,
                                               Vector2& normalized_key)
  {
    Scalar theta = rotation.norm();
    Vector3 camera_point;
    if (theta == Scalar(0))
    {
      camera_point = world_point + translation;
    }
    else
    {
      Vector3 normalized_rotation = rotation / theta;
      camera_point =
        std::cos(theta) * world_point +
        std::sin(theta) * normalized_rotation.cross(world_point) +
        (1 - std::cos(theta)) * world_point.dot(normalized_rotation) *
        normalized_rotation + translation;
    }
    camera_point /= camera_point[2];
    normalized_key = camera_point.segment(0, 2);
  }

  inline static void GetRadialDistortDelta(Scalar k1,
                                           Scalar k2,
                                           Scalar k3,
                                           const Vector2& normalized_key,
                                           Vector2& delta)
  {
    RadialDistortor<Scalar>()(k1, k2, k3, normalized_key[0], normalized_key[1],
                              delta[0], delta[1]);
  }

  inline static void GetDecenteringDistortDelta(Scalar d1,
                                                Scalar d2,
                                                const Vector2& normalized_key,
                                                Vector2& delta)
  {
    DecenteringDistortor<Scalar>()(d1, d2, normalized_key[0], normalized_key[1],
                                   delta[0], delta[1]);
  }

  inline static void NormalizedKeyToImageKey(Scalar focal_length,
                                             Scalar skew,
                                             Scalar principal_x,
                                             Scalar principal_y,
                                             Scalar pixel_ratio,
                                             const Vector2& normalized_key,
                                             Vector2& image_key)
  {
    image_key[0] = focal_length * normalized_key[0] +
                   skew * normalized_key[1] +
                   principal_x;
    image_key[1] = focal_length * pixel_ratio * normalized_key[1] +
                   principal_y;
  }

private:
  Index number_of_images_;
  Index number_of_points_;
  Index number_of_keys_;
  Index number_of_cameras_;
  FeatureMapContainer feature_maps_;
  ImageCameraMap image_camera_map_;

  IntrinsicComputationsMask intrinsic_computations_mask_;

  PointConstraintContainer point_constraints_;
  ImageConstraintContainer image_constraints_;
  CameraConstraintContainer camera_constraints_;

  FixMask fix_mask_;
  PointContainer fix_points_;
  ImageContainer fix_images_;
  CameraContainer fix_cameras_;
};

}
}
}

#endif
