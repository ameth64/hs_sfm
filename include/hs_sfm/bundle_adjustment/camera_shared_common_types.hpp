#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_COMMON_TYPES_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_COMMON_TYPES_HPP_

#include <bitset>
#include <vector>

namespace hs
{
namespace sfm
{
namespace ba
{

enum
{
  COMPUTE_RADIAL_DISTORTION = 0,
  COMPUTE_DECENTERING_DISTORTION,
  COMPUTE_INTRINSIC_PARAMS,
  NUMBER_OF_INTRINSIC_COMPUTATIONS
};
typedef std::bitset<NUMBER_OF_INTRINSIC_COMPUTATIONS>
        IntrinsicComputationsMask;

enum
{
  POINT_CONSTRAIN_X = 0,
  POINT_CONSTRAIN_Y,
  POINT_CONSTRAIN_Z,
  NUMBER_OF_POINT_CONSTRAINTS
};
typedef std::bitset<NUMBER_OF_POINT_CONSTRAINTS>
        PointConstraintsMask;
struct PointConstraint
{
  size_t point_id;
  PointConstraintsMask mask;
};
typedef std::vector<PointConstraint> PointConstraintContainer;

enum
{
  IMAGE_CONSTRAIN_ROTATION = 0,
  IMAGE_CONSTRAIN_POSITION_X,
  IMAGE_CONSTRAIN_POSITION_Y,
  IMAGE_CONSTRAIN_POSITION_Z,
  NUMBER_OF_IMAGE_CONSTRAINTS
};
typedef std::bitset<NUMBER_OF_IMAGE_CONSTRAINTS>
        ImageConstraintsMask;
struct ImageConstraint
{
  size_t image_id;
  ImageConstraintsMask mask;
};
typedef std::vector<ImageConstraint> ImageConstraintContainer;

enum
{
  RADIAL_CONSTRAIN_K1 = 0,
  RADIAL_CONSTRAIN_K2,
  RADIAL_CONSTRAIN_K3,
  NUMBER_OF_RADIAL_CONSTRAINTS
};
typedef std::bitset<NUMBER_OF_RADIAL_CONSTRAINTS>
        RadialConstraintsMask;

enum
{
  DECENTERING_CONSTRAIN_D1 = 0,
  DECENTERING_CONSTRAIN_D2,
  NUMBER_OF_DECENTERING_CONSTRAINTS
};
typedef std::bitset<NUMBER_OF_DECENTERING_CONSTRAINTS>
        DecenteringConstraintsMask;

enum
{
  INTRINSIC_CONSTRAIN_FOCAL_LENGTH = 0,
  INTRINSIC_CONSTRAIN_SKEW,
  INTRINSIC_CONSTRAIN_PRINCIPAL_X,
  INTRINSIC_CONSTRAIN_PRINCIPAL_Y,
  INTRINSIC_CONSTRAIN_PIXEL_RATIO,
  NUMBER_OF_INTRINSIC_CONSTRAINTS
};
typedef std::bitset<NUMBER_OF_INTRINSIC_CONSTRAINTS>
        IntrinsicConstraintsMask;

struct CameraConstraint
{
  size_t camera_id;
  RadialConstraintsMask radial_mask;
  DecenteringConstraintsMask decentering_mask;
  IntrinsicConstraintsMask intrinsic_mask;
};
typedef std::vector<CameraConstraint> CameraConstraintContainer;

}
}
}

#endif
