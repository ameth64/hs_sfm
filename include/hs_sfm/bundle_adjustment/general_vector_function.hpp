#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_GENERAL_VECTOR_FUNCTION_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_GENERAL_VECTOR_FUNCTION_HPP_

#include <vector>
#include <utility>
#include <bitset>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/sfm_utility/radial_distortor.hpp"
#include "hs_sfm/sfm_utility/decentering_distortor.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class GeneralVectorFunction
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

  enum
  {
    COMPUTE_RADIAL_DISTORTION = 0,
    COMPUTE_DECENTERING_DISTORTION,
    COMPUTE_INTRINSIC_PARAMS,
    NUMBER_OF_INTRINSIC_COMPUTATION
  };

  enum
  {
    CONSTRAIN_RADIAL_K1 = 0,
    CONSTRAIN_RADIAL_K2,
    CONSTRAIN_RADIAL_K3,
    CONSTRAIN_DECENTERING_D1,
    CONSTRAIN_DECENTERING_D2,
    CONSTRAIN_FOCAL_LENGTH,
    CONSTRAIN_SKEW,
    CONSTRAIN_PRINCIPAL_X,
    CONSTRAIN_PRINCIPAL_Y,
    CONSTRAIN_PIXEL_RATIO,
    NUMBER_OF_INTRINSIC_CONSTRAINTS
  }

  enum
  {
    CONSTRAIN_POINTS = 0,
    NUMBER_OF_STRUCTURE_CONSTRAINTS;
  }

  typedef std::bitset<NUMBER_OF_INTRINSIC_COMPUTATIONS>
          IntrinsicComputationsMask;
  typedef std::bitset<NUMBER_OF_INTRINSIC_CONSTRAINTS>
          IntrinsicConstraintsMask;
  typedef std::bitset<NUMBER_OF_STRUCTURE_CONSTRAINTS>
          StructureConstraintsMask;

  typedef std::pair<Index, Index> FeatureMap;
  typedef std::vector<FeatureMap> FeatureMapContainer;

  typedef std::vector<Index> ImageCameraMap;

  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;

public:
  GeneralVectorFunction()
    : number_of_images_(0),
      number_of_points_(0),
      number_of_keys_(0),
      number_of_constrained_points_(0),
      number_of_cameras_(0),
      feature_maps_(),
      image_camera_map_(),
      intrinsic_computations_mask_(),
      intrinsic_constraints_mask_(),
      structure_constraints_mask_() {}

  inline Index number_of_cameras() const
  {
    return number_of_cameras_;
  }
  inline void set_number_of_cameras(Index number_of_cameras)
  {
    number_of_cameras_ = number_of_cameras_;
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

  inline Index number_of_constrained_points() const
  {
    return number_of_constrained_points_;
  }
  inline void set_number_of_constrained_points(
    Index number_of_constrained_points)
  {
    number_of_constrained_points_ = number_of_constrained_points;
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

  inline const IntrinsicConstraintsMask& intrinsic_constraints_mask() const
  {
    return intrinsic_constraints_mask_;
  }
  inline IntrinsicConstraintsMask& intrinsic_constraints_mask()
  {
    return intrinsic_constraints_mask_;
  }

  inline const StructureConstraintsMask& structure_constraints_mask() const
  {
    return structure_constraints_mask_;
  }
  inline StructureConstraintsMask& structure_constraints_mask()
  {
    return structure_constraints_mask_;
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

  inline Index GetYConstrainedPointsSize() const
  {
    if (structure_constraints_mask_[CONSTRAIN_POINTS])
    {
      return number_of_constrained_points_ * params_per_point_;
    }
    else
    {
      return 0;
    }
  }

  inline Index GetYConstrainedIntrinsicParamsSize() const
  {
    Index intrinsic_params_per_camera =
      Index(intrinsic_constraints_mask_.count());
    return intrinsic_params_per_camera * number_of_cameras_;
  }

  inline Index GetXSize() const
  {
    return (GetExtrinsicParamsSize() +
            GetPointParamsSize() +
            GetIntrinsicParamsSize());
  }

  inline Index GetYSize() const
  {
    return (GetYKeysSize() +
            GetYConstrainedPointsSize() +
            GetYConstrainedIntrinsicParamsSize());
  }

  Err operator() (const XVector& x, YVector& y) const
  {
    Index x_size = x.rows();
    if (x_size != GetXSize()) return -1;
    if (number_of_keys_ != Index(feature_maps_.size())) return -1;
    if (number_of_images_ != Index(image_camera_map_.size())) return -1;
    if ((intrinsic_constraints_mask_[CONSTRAIN_FOCAL_LENGTH] ||
         intrinsic_constraints_mask_[CONSTRAIN_SKEW] ||
         intrinsic_constraints_mask_[CONSTRAIN_PRINCIPAL_X] ||
         intrinsic_constraints_mask_[CONSTRAIN_PRINCIPAL_Y] ||
         intrinsic_constraints_mask_[CONSTRAIN_PIXEL_RATIO]) &&
        !intrinsic_computations_mask_[COMPUTE_INTRINSIC_PARAMS]) return -1;
    if ((intrinsic_constraints_mask_[CONSTRAIN_RADIAL_K1] ||
         intrinsic_constraints_mask_[CONSTRAIN_RADIAL_K2] ||
         intrinsic_constraints_mask_[CONSTRAIN_RADIAL_K3]) &&
        !intrinsic_computations_mask_[COMPUTE_RADIAL_DISTORTION]) return -1;
    if ((intrinsic_constraints_mask_[CONSTRAIN_DECENTERING_D1] ||
         intrinsic_constraints_mask_[CONSTRAIN_DECENTERING_D2]) &&
        !intrinsic_computations_mask_[COMPUTE_DECENTERING_DISTORTION]) return -1;
    Index y_size = GetYSize();
    y.resize(y_size);

    //keys segment
    Index point_params_size = GetPointParamsSize();
    Index extrinsic_params_size = GetExtrinsicParamsSize();
    Index intrinsic_params_size_per_camera =
      GetIntrinsicParamsSizePerCamera();
    for (Index i = 0; i < number_of_keys_; i++)
    {
      Index image_id = feature_maps_[i].first;
      Index point_id = feature_maps_[i].second;

      Vector3 point = x.segment(point_id * params_per_point_,
                                params_per_point_);
      Vector3 rotation =
        x.segment(point_params_size + image_id * params_per_camera, 3);
      Vector3 translation =
        x.segment(point_params_size + image_id * params_per_camera, + 3, 3);

      Vector2 normalized_key;
      WorldPointToNormalizedKey(point,
                                rotation,
                                translation,
                                normalized_key);

      Index camera_id = image_camera_map_[image_id];
      Index offset = point_params_size + extrinsic_params_size +
                     camera_id * intrinsic_params_size_per_camera;
      if (intrinsic_computations_mask_[COMPUTE_RADIAL_DISTORTION])
      {
        Scalar k1 = x[offset + 0];
        Scalar k2 = x[offset + 1];
        Scalar k3 = x[offset + 2];
        RadialDistortNormalizedKey(k1, k2, k3, normalized_key);
      }
      if (intrinsic_computations_mask_[COMPUTE_DECENTERING_DISTORTION])
      {
        offset +=
          intrinsic_computations_mask_[COMPUTE_RADIAL_DISTORTION] ? 3 : 0;
        Scalar d1 = x[offset + 0];
        Scalar d2 = x[offset + 1];
        DecenteringDistortNormalizedKey(d1, d2, normalized_key);
      }
      Vector2 image_key;
      if (intrinsic_computations_mask_[COMPUTE_INTRINSIC_PARAMS])
      {
        offset +=
          intrinsic_computations_mask_[COMPUTE_DECENTERING_DISTORTION] ? 2 : 0;
        Scalar focal_length = x[offset + 0];
        Scalar skew = x[offset + 1];
        Scalar principal_x = x[offset + 2];
        Scalar principal_y = x[offset + 3];
        Scalar pixel_ratio = x[offset + 4];
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

      y.segment(i * params_per_key_, params_per_key_) = image_key;
    }

    Index y_offset = GetYKeysSize();
    if (intrinsic_constraints_mask_.any())
    {
      for (Index i = 0; i < number_of_cameras_; i++)
      {
        Index x_intrinsic_offset = x_offset +
                                   i * intrinsic_params_size_per_camera;
        Index y_intrinsic_offset =
          y_offset + i * Index(intrinsic_constraints_mask_.count());
        if (intrinsic_constraints_mask_[CONSTRAIN_RADIAL_K1])
        {
          y[y_intrinsic_offset] = x[x_intrinsic_offset + 0];
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[CONSTRAIN_RADIAL_K2])
        {
          y[y_intrinsic_offset] = x[x_intrinsic_offset + 1];
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[CONSTRAIN_RADIAL_K3])
        {
          y[y_intrinsic_offset] = x[x_intrinsic_offset + 2];
          y_intrinsic_offset++;
        }
        x_intrinsic_offset +=
          intrinsic_computations_mask_[COMPUTE_RADIAL_DISTORTION] ? 3 : 1;
        if (intrinsic_constraints_mask_[CONSTRAIN_DECENTERING_D1])
        {
          y[y_intrinsic_offset] = x[x_intrinsic_offset + 0];
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[CONSTRAIN_DECENTERING_D2])
        {
          y[y_intrinsic_offset] = x[x_intrinsic_offset + 1];
          y_intrinsic_offset++;
        }
        x_intrinsic_offset +=
          intrinsic_computations_mask_[COMPUTE_DECENTERING_DISTORTION] ? 2 : 1;
        if (intrinsic_constraints_mask_[CONSTRAIN_FOCAL_LENGTH])
        {
          y[y_intrinsic_offset] = x[x_intrinsic_offset + 0];
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[CONSTRAIN_SKEW])
        {
          y[y_intrinsic_offset] = x[x_intrinsic_offset + 1];
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[CONSTRAIN_PRINCIPAL_X])
        {
          y[y_intrinsic_offset] = x[x_intrinsic_offset + 2];
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[CONSTRAIN_PRINCIPAL_Y])
        {
          y[y_intrinsic_offset] = x[x_intrinsic_offset + 3];
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[CONSTRAIN_PIXEL_RATIO])
        {
          y[y_intrinsic_offset] = x[x_intrinsic_offset + 4];
          y_intrinsic_offset++;
        }
      }// for (Index i = 0; i < number_of_cameras_; i++)
      y_offset += GetYConstrainedIntrinsicParamsSize();
    }// if (intrinsic_constraints_mask_.any())

    if (structure_constraints_mask_.any())
    {
      if (structure_constraints_mask_[CONSTRAIN_POINTS])
      {
        Index x_offset = GetExtrinsicParamsSize() +
                         (number_of_points_ - number_of_constrained_points_) *
                         params_per_point_;
        for (Index i = 0; i < number_of_constrained_points_; i++)
        {
          Index y_point_offset = y_offset + i * params_per_point_;
          Index x_point_offset = x_offset + i * params_per_point_;
          y.segment(y_point_offset, params_per_point_) =
            x.segment(x_point_offset, params_per_point_);
        }
      }// if (structure_constraints_mask_[CONSTRAIN_POINTS])
    }// if (structure_constraints_mask_.any())

    return 0;
  }

  inline static void WorldPointToNormalizedKey(const Vecotr3& world_point,
                                               const Vector3& rotation,
                                               const Vector3& translation,
                                               Vector2& normalized_key)
  {
    Scalar theta = rotation.norm();
    vector3 camera_point;
    if (theta == Scalar(0))
    {
      camera_point = point + translation;
    }
    else
    {
      Vector3 normalized_rotation = rotation / theta;
      camera_point =
        std::cos(theta) * point +
        std::sin(theta) * normalized_rotation.cross(point) +
        (1 - std::cos(theta)) * point.dot(normalized_rotation) +
        normalized_rotation + translation;
    }
    camera_point /= camera_point[2];
    normalized_key = camera_point.segment(0, 2);
  }

  inline static void RadialDistortNormalizedKey(Scalar k1,
                                                Scalar k2,
                                                Scalar k3,
                                                Vector2& normalized_key)
  {
    Vector2 delta_key;
    RadialDistortor()(k1, k2, k3, normalized_key[0], normalized_key[1],
                      delta_key[0], delta_key[1]);
    normalized_key += delta_key;
  }

  inline static void DecenteringDistortNormalizedKey(Scalar d1,
                                                     Scalar d2,
                                                     Vector2& normalized_key)
  {
    Vector2 delta_key;
    DecenteringDistortor()(d1, d2, normalized_key[0], normalized_key[1],
                           delta_key[0], delta_key[1]);
    normalized_key += delta_key;
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
  Index number_of_constrained_points_;
  Index number_of_cameras_;
  FeatureMapContainer feature_maps_;
  ImageCameraMap image_camera_map_;

  IntrinsicComputationsMask intrinsic_computations_mask_;
  IntrinsicConstraintsMask intrinsic_constraints_mask_;
  StructureConstraintsMask structure_constraints_mask_;
};

}
}
}

#endif
