#ifndef _HS_SFM_SYNTHETIC_KEYSET_GENERATOR_HPP_
#define _HS_SFM_SYNTHETIC_KEYSET_GENERATOR_HPP_

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
class KeysetGenerator
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
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;

  typedef int Err;

  KeysetGenerator(ImageDimension image_width, ImageDimension image_height)
    : image_width_(image_width), image_height_(image_height) {}

  Err operator()(const IntrinContainer& intrinsic_params_set,
                 const ExtrinContainer& extrinsic_params_set,
                 const Point3DContainer& point,
                 KeysetContainer& keysets,
                 TrackContainer& tracks,
                 CameraViewContainer& camera_views) const
  {
    size_t number_of_cameras = intrinsic_params_set.size();
    if (number_of_cameras != extrinsic_params_set.size())
    {
    return -1;
    }

    Point2DSetContainer point2dsets(number_of_cameras);
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
            std::make_pair(j, point2dsets[j].size()));
          camera_views[j].push_back(
            std::make_pair(i, point2dsets[j].size()));
          point2dsets[j].push_back(key_homogeneous.segment(0, 2));
        }
      }
    }
    keysets.clear();
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      keysets.push_back(Keyset(point2dsets[i]));
    }

    return 0;
  }

  ImageDimension image_width() const
  {
    return image_width_;
  }

  ImageDimension image_height() const
  {
    return image_height_;
  }

  private:
  ImageDimension image_width_;
  ImageDimension image_height_;
};

}
}
}

#endif
