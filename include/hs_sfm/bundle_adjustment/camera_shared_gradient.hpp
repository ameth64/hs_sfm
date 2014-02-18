#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_GRADIENT_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_GRADIENT_HPP_

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
              Index camera_params_size)
  {
    number_of_points_ = number_of_points;
    number_of_images_ = number_of_images;
    number_of_cameras_ = number_of_cameras;
    camera_params_size_ = camera_params_size;

    point_segments_.resize(number_of_points_, PointSegment::Zero());
    image_segments_.resize(number_of_images_, ImageSegment::Zero());
    camera_segments_.resize(number_of_cameras_,
                            CameraSegment::Zero(camera_params_size_));
  }

  Scalar operator[] (Index i) const
  {
    Index points_end = Index(number_of_points_) * params_per_point_;
    Index images_end =
      points_end + Index(number_of_images_) * extrinsic_params_per_image_;
    Index cameras_end =
      images_end + Index(number_of_cameras_) * camera_params_size_;
    if (i < 0 || i >= cameras_end)
    {
      return std::numeric_limits<Scalar>::signaling_NaN();
    }
    else if(i < points_end)
    {
      Index point_id = i / params_per_point_;
      Index param_id = i % params_per_point_;
      return point_segments_[point_id][param_id];
    }
    else if (i < images_end)
    {
      Index image_id = (i - points_end) / extrinsic_params_per_image_;
      Index param_id = (i - points_end) % extrinsic_params_per_image_;
      return image_segments_[image_id][param_id];
    }
    else
    {
      Index camera_id = (i - images_end) / camera_params_size_;
      Index param_id = (i - images_end) % camera_params_size_;
      return camera_segments_[camera_id][param_id];
    }
  }

  const PointSegment& GetPointSegment(Index point_id) const
  {
    return point_segments_[point_id];
  }
  PointSegment& GetPointSegment(Index point_id)
  {
    return point_segments_[point_id];
  }

  const ImageSegment& GetImageSegment(Index image_id) const
  {
    return image_segments_[image_id];
  }
  ImageSegment& GetImageSegment(Index image_id)
  {
    return image_segments_[image_id];
  }

  const CameraSegment& GetCameraSegment(Index camera_id) const
  {
    return camera_segments_[camera_id];
  }
  CameraSegment& GetCameraSegment(Index camera_id)
  {
    return camera_segments_[camera_id];
  }

  Index GetXSize() const
  {
    return Index(number_of_points_) * params_per_point_ +
           Index(number_of_images_) * extrinsic_params_per_image_ +
           Index(number_of_cameras_) * camera_params_size_;
  }

private:
  size_t number_of_points_;
  size_t number_of_images_;
  size_t number_of_cameras_;
  Index camera_params_size_;
  PointSegmentContainer point_segments_;
  ImageSegmentContainer image_segments_;
  CameraSegmentContainer camera_segments_;
};

}
}
}

#endif
