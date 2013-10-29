#ifndef _HS_SFM_SYNTHETIC_RELATIVE_GENERATOR_HPP_
#define _HS_SFM_SYNTHETIC_RELATIVE_GENERATOR_HPP_

namespace hs
{
namespace sfm
{
namespace synthetic
{

template <typename _Scalar>
class RelativeGenerator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef typename ExtrinsicParams::Rotation rotation;
  typedef typename ExtrinsicParams::Position Position;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;
  typedef ImageKeys<Scalar> ImageKeyset;
  typedef EIGEN_STD_VECTOR(ImageKeyset) ImageKeysetContainer;

  typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;
  typedef EIGEN_VECTOR(Scalar, 3) Translate;

private:
  class CompareFunctor
  {
  public:
    CompareFunctor(const ExtrinsicParamsContainer& cameras,
                   const ExtrinsicParams& fix_camera)
      : cameras_(cameras),
        fix_camera_(fix_camera) {}

    bool operator() (size_t id_1, size_t id_2) const
    {
      Scalar distance1 = (cameras_[id_1].position() -
                          fix_camera_.position()).norm();
      Scalar distance2 = (cameras_[id_2].position() -
                          fix_camera_.position()).norm();
      return (distance1 < distance2);
    }

    const ExtrinsicParamsContainer& cameras_;
    const ExtrinsicParams& fix_camera_;
  };

public:
  Err operator() (const TrackContainer& tracks,
                  const CameraViewContainer& camera_views,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  Point3DContainer& points,
                  RMatrix& rotation_similar,
                  Translate& translate_similar,
                  Scalar& scale_similar,
                  size_t& camera_id_identity,
                  size_t& camera_id_relative) const
  {
    if (ChooseCameraPair(extrinsic_params_set, camera_views,
                         camera_id_identity,
                         camera_id_relative) != 0)
    {
      return -1;
    }

    if (ComputeSimilarTransform(extrinsic_params_set,
                                camera_id_identity,
                                camera_id_relative,
                                rotation_similar,
                                translate_similar,
                                scale_similar) != 0)
    {
      return -1;
    }

    if (SimilarTransformCameras(rotation_similar,
                                translate_similar,
                                scale_similar,
                                extrinsic_params_set) != 0)
    {
      return -1;
    }

    if (SimilarTransformPoints(rotation_similar,
                               translate_similar,
                               scale_similar,
                               points) != 0)
    {
      return -1;
    }

    return 0;
  }

private:
  Err ChooseCameraPair(const ExtrinsicParamsContainer& extrinsic_params_set,
                       const CameraViewContainer& camera_views,
                       size_t& camera_id_identity,
                       size_t& camera_id_relative) const
  {
    size_t number_of_cameras = camera_views.size();
    if (number_of_cameras != extrinsic_params_set.size() ||
        number_of_cameras < 2)
    {
      return -1;
    }

    //选取包含特征点最多的相机作为固定的相机
    size_t max_number_of_camera_views = camera_views[0].size();
    camera_id_identity = 0;
    for (size_t i = 1; i < number_of_cameras; i++)
    {
      if (camera_views[i].size() > max_number_of_camera_views)
      {
        max_number_of_camera_views = camera_views[i].size();
        camera_id_identity = i;
      }
    }

    //选取距离固定相机第四近的相机作为相对相机
    CompareFunctor compare_functor(extrinsic_params_set,
                                   extrinsic_params_set[camera_id_identity]);
    std::vector<size_t> indices(number_of_cameras);
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(), compare_functor);

    if (number_of_cameras < 4)
    {
      camera_id_relative = indices[number_of_cameras - 1];
    }
    else
    {
      camera_id_relative = indices[3];
    }

    return 0;
  }

  Err ComputeSimilarTransform(
    const ExtrinsicParamsContainer& extrinsic_params_set,
    size_t camera_id_identity,
    size_t camera_id_relative,
    RMatrix& rotation_similar,
    Translate& translate_similar,
    Scalar& scale_similar) const
  {
    RMatrix rotation_identity =
      extrinsic_params_set[camera_id_identity].rotation();

    Position position_identity =
      extrinsic_params_set[camera_id_identity].position();
    Position position_relative =
      extrinsic_params_set[camera_id_relative].position();

    rotation_similar = rotation_identity.transpose();
    translate_similar = position_identity;
    scale_similar = (position_relative - position_identity).norm();

    return 0;
  }

  Err SimilarTransformCameras(
    const RMatrix& rotation_similar,
    const Translate& translate_similar,
    Scalar scale_similar,
    ExtrinsicParamsContainer& extrinsic_params_set) const
  {
    size_t number_of_cameras = extrinsic_params_set.size();
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      RMatrix rotation_absolute = extrinsic_params_set[i].rotation();
      Position position_absolute = extrinsic_params_set[i].position();

      RMatrix rotation_relative = rotation_absolute * rotation_similar;
      Position position_relative =
        rotation_similar.transpose() *
        (position_absolute - translate_similar) / scale_similar;

      extrinsic_params_set[i].rotation() = rotation_relative;
      extrinsic_params_set[i].position() = position_relative;
    }
    return 0;
  }

  Err SimilarTransformPoints(
    const RMatrix& rotation_similar,
    const Translate& translate_similar,
    Scalar scale_similar,
    Point3DContainer& points) const
  {
    size_t number_of_points = points.size();
    for (size_t i = 0; i < number_of_points; i++)
    {
      points[i] =
        rotation_similar.transpose() * (points[i] -
                                        translate_similar) / scale_similar;
    }
    return 0;
  }

};

}
}
}

#endif
