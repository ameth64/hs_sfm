#ifndef _BA_GCP_CONSTRAINED_RELATIVE_DATA_GENERATOR_HPP_
#define _BA_GCP_CONSTRAINED_RELATIVE_DATA_GENERATOR_HPP_

#include <vector>
#include <algorithm>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_math/geometry/rotation.hpp"

#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BAGCPConstrainedRelativeDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef BAGCPConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::FeatureMap FeatureMap;
  typedef typename VectorFunction::FeatureMapContainer FeatureMapContainer;
  typedef hs::math::geometry::Rotation3D<Scalar> Rotation;
  typedef EIGEN_VECTOR(Scalar, 3) Translate;

private:
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;
  struct CompareFunctor
  {
    bool operator() (size_t id_1, size_t id_2) const
    {
      Scalar distance1 = (points[id_1] - fix_point).norm();
      Scalar distance2 = (points[id_2] - fix_point).norm();
      return (distance1 < distance2);
    }

    Point fix_point;
    PointContainer points;
  };

public:
  Err operator()(const VectorFunction& vector_function,
                 const XVector& absolute_x,
                 XVector& relative_x,
                 Scalar& scale,
                 Rotation& rotation,
                 Translate& translate) const
  {
    Index camera_id_identity;
    Index camera_id_relative;
    if (ChooseCameraPair(vector_function, absolute_x,
                         camera_id_identity, camera_id_relative) != 0)
    {
      return -1;
    }

    if (SimilarTransformXVector(vector_function, absolute_x,
                                camera_id_identity, camera_id_relative,
                                relative_x, scale, rotation, translate) != 0)
    {
      return -1;
    }

    return 0;
  }

private:
  Err ChooseCameraPair(const VectorFunction& vector_function,
                       const XVector& absolute_x,
                       Index& camera_id_identity,
                       Index& camera_id_relative) const
  {
    //选取包含特征点最多的相机作为固定的相机
    size_t number_of_cameras = size_t(vector_function.number_of_cameras());
    std::vector<size_t> number_of_camera_features(number_of_cameras);
    const FeatureMapContainer& feature_maps = vector_function.feature_maps();
    for (size_t i = 0; i < feature_maps.size(); i++)
    {
      number_of_camera_features[feature_maps[i].first]++;
    }
    size_t max_number_of_camera_features = 0;
    Index max_features_camera_id = -1;
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      if (max_number_of_camera_features < number_of_camera_features[i])
      {
        max_number_of_camera_features = number_of_camera_features[i];
        max_features_camera_id = i;
      }
    }

    camera_id_identity = max_features_camera_id;

    //选取距离固定相机第四近的相机作为相对相机
    CompareFunctor compare_functor;

    compare_functor.points.resize(number_of_cameras);
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      Point rotation_vector =
        absolute_x.segment(i * VectorFunction::params_per_camera_, 3);
      Rotation rotation(rotation_vector);
      Point camera_t = 
        absolute_x.segment(i * VectorFunction::params_per_camera_ + 3, 3);
      Point camera_c = -(rotation.Inverse() * camera_t);
      compare_functor.points[i] = camera_c;
    }
    compare_functor.fix_point = compare_functor.points[camera_id_identity];
    std::vector<size_t> indices(number_of_cameras);
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(), compare_functor);

    camera_id_relative = Index(indices[3]);

    return 0;
  }

  Err SimilarTransformXVector(const VectorFunction& vector_function,
                              const XVector& absolute_x,
                              Index camera_id_identity,
                              Index camera_id_relative,
                              XVector& relative_x,
                              Scalar& scale,
                              Rotation& rotation,
                              Translate& translate) const
  {
    Point identity_rotation_vector =
      absolute_x.segment(
        camera_id_identity * VectorFunction::params_per_camera_, 3);
    Rotation identity_rotation(identity_rotation_vector);
    Point identity_t =
      absolute_x.segment(
        camera_id_identity * VectorFunction::params_per_camera_ + 3, 3);
    Point identity_c = -(identity_rotation.Inverse() * identity_t);

    Point relative_rotation_vector =
      absolute_x.segment(
        camera_id_relative * VectorFunction::params_per_camera_, 3);
    Rotation relative_rotation(relative_rotation_vector);
    Point relative_t =
      absolute_x.segment(
        camera_id_relative * VectorFunction::params_per_camera_ + 3, 3);
    Point relative_c = -(relative_rotation.Inverse() * relative_t);

    rotation = identity_rotation.Inverse();
    translate = -(identity_rotation.Inverse() * identity_t);
    scale = (relative_c - identity_c).norm();

    Index x_size = vector_function.GetXSize();
    relative_x.resize(x_size);
    Index number_of_cameras = vector_function.number_of_cameras();
    for (Index i = 0; i < number_of_cameras; i++)
    {
      Point absolute_camera_rotation_vector =
        absolute_x.segment(
          i * VectorFunction::params_per_camera_, 3);
      Point absolute_camera_t =
        absolute_x.segment(
          i * VectorFunction::params_per_camera_ + 3, 3);
      Rotation absolute_camera_rotation(absolute_camera_rotation_vector);
      Rotation relative_camera_rotation = absolute_camera_rotation * rotation;
      Point relative_camera_t = absolute_camera_t +
                                absolute_camera_rotation * translate;
      relative_camera_t /= scale;
      relative_x[i * VectorFunction::params_per_camera_ + 0] =
        relative_camera_rotation[0];
      relative_x[i * VectorFunction::params_per_camera_ + 1] =
        relative_camera_rotation[1];
      relative_x[i * VectorFunction::params_per_camera_ + 2] =
        relative_camera_rotation[2];
      relative_x.segment(i * VectorFunction::params_per_camera_ + 3, 3) =
        relative_camera_t;
    }
    
    Index number_of_points = vector_function.number_of_points();
    Index point_start_id = vector_function.GetCameraParamsSize();
    for (Index i = 0; i < number_of_points; i++)
    {
      Point absolute_point =
        absolute_x.segment(point_start_id + 
                           i * VectorFunction::params_per_point_,
                           VectorFunction::params_per_point_);
      Point relative_point = (absolute_point - translate) / scale;
      relative_point = rotation.Inverse() * relative_point;
      relative_x.segment(point_start_id + i * VectorFunction::params_per_point_,
                         VectorFunction::params_per_point_) = relative_point;
    }

    return 0;
  }
};

}
}
}

#endif