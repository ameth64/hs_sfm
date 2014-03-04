#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_GRADIENT_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_GRADIENT_HPP_

#include <cassert>

#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedGradient
{
public:
  typedef _Scalar Scalar;

  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;

  static const Index params_per_key_ = VectorFunction::params_per_key_;
  static const Index params_per_point_ = VectorFunction::params_per_point_;
  static const Index extrinsic_params_per_image_ =
    VectorFunction::extrinsic_params_per_image_;

  typedef EIGEN_VECTOR(Scalar, params_per_point_) PointSegment;
  typedef EIGEN_VECTOR(Scalar, extrinsic_params_per_image_) ImageSegment;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) CameraSegment;
  typedef EIGEN_STD_VECTOR(PointSegment) PointSegmentContainer;
  typedef EIGEN_STD_VECTOR(ImageSegment) ImageSegmentContainer;
  typedef EIGEN_STD_VECTOR(CameraSegment) CameraSegmentContainer;
public:
  void Resize(size_t number_of_points,
              size_t number_of_images,
              size_t number_of_cameras,
              Index camera_params_size,
              const FixMask& fix_mask)
  {
    number_of_points_ = number_of_points;
    number_of_images_ = number_of_images;
    number_of_cameras_ = number_of_cameras;
    camera_params_size_ = camera_params_size;

    fix_mask_ = fix_mask;

    if (!fix_mask_[FIX_POINTS])
      point_segments_.resize(number_of_points_, PointSegment::Zero());
    if (!fix_mask_[FIX_IMAGES])
      image_segments_.resize(number_of_images_, ImageSegment::Zero());
    if (!fix_mask_[FIX_CAMERAS])
      camera_segments_.resize(number_of_cameras_,
                              CameraSegment::Zero(camera_params_size_));
  }

  Scalar operator[] (Index i) const
  {
    Index x_size = GetXSize();
    Index points_end = fix_mask_[FIX_POINTS] ?
      0 : (Index(number_of_points_) * params_per_point_);
    Index images_end = points_end + (fix_mask_[FIX_IMAGES] ?
       0 : Index(number_of_images_) * extrinsic_params_per_image_);
    Index cameras_end = images_end + (fix_mask_[FIX_CAMERAS] ?
       0 : Index(number_of_cameras_) * camera_params_size_);
    if (i < 0 || i >= x_size)
    {
      return std::numeric_limits<Scalar>::signaling_NaN();
    }
    else if(i < points_end && !fix_mask_[FIX_POINTS])
    {
      Index point_id = i / params_per_point_;
      Index param_id = i % params_per_point_;
      return point_segments_[point_id][param_id];
    }
    else if (i < images_end && !fix_mask_[FIX_IMAGES])
    {
      Index image_id = (i - points_end) / extrinsic_params_per_image_;
      Index param_id = (i - points_end) % extrinsic_params_per_image_;
      return image_segments_[image_id][param_id];
    }
    else if (i < cameras_end && !fix_mask_[FIX_CAMERAS])
    {
      Index camera_id = (i - images_end) / camera_params_size_;
      Index param_id = (i - images_end) % camera_params_size_;
      return camera_segments_[camera_id][param_id];
    }
    else
    {
      return std::numeric_limits<Scalar>::signaling_NaN();
    }
  }

  const PointSegment& GetPointSegment(Index point_id) const
  {
    assert(!fix_mask_[FIX_POINTS]);
    return point_segments_[point_id];
  }
  PointSegment& GetPointSegment(Index point_id)
  {
    assert(!fix_mask_[FIX_POINTS]);
    return point_segments_[point_id];
  }

  const ImageSegment& GetImageSegment(Index image_id) const
  {
    assert(!fix_mask_[FIX_IMAGES]);
    return image_segments_[image_id];
  }
  ImageSegment& GetImageSegment(Index image_id)
  {
    assert(!fix_mask_[FIX_IMAGES]);
    return image_segments_[image_id];
  }

  const CameraSegment& GetCameraSegment(Index camera_id) const
  {
    assert(!fix_mask_[FIX_CAMERAS]);
    return camera_segments_[camera_id];
  }
  CameraSegment& GetCameraSegment(Index camera_id)
  {
    assert(!fix_mask_[FIX_CAMERAS]);
    return camera_segments_[camera_id];
  }

  Index GetXSize() const
  {
    return (fix_mask_[FIX_POINTS] ?
           0 : (Index(number_of_points_) * VectorFunction::params_per_point_)) +
           (fix_mask_[FIX_IMAGES] ?
           0 : (Index(number_of_images_) *
                VectorFunction::extrinsic_params_per_image_)) +
           (fix_mask_[FIX_CAMERAS] ?
           0 : (Index(number_of_cameras_) * camera_params_size_));
  }

private:
  size_t number_of_points_;
  size_t number_of_images_;
  size_t number_of_cameras_;
  Index camera_params_size_;
  PointSegmentContainer point_segments_;
  ImageSegmentContainer image_segments_;
  CameraSegmentContainer camera_segments_;
  FixMask fix_mask_;
};

}
}
}

#endif
