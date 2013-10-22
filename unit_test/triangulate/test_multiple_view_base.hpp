#ifndef _UNIT_TEST_TEST_MULTIPLE_VIEW_BASE_HPP_
#define _UNIT_TEST_TEST_MULTIPLE_VIEW_BASE_HPP_

#include <gtest/gtest.h>

#include "hs_sfm/sfm_utility/synthetic_scene_generator.hpp"
#include "hs_sfm/triangulate/multiple_view_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace triangulate
{

template <typename _Scalar, typename _ImageDimension>
class MultipleViewSyntheicDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;

  typedef int Err;

  typedef hs::sfm::SceneGenerator<Scalar, ImageDimension> SceneGenerator;
  typedef typename SceneGenerator::IntrinsicParams IntrinsicParams;
  typedef typename SceneGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename SceneGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename SceneGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename SceneGenerator::Image Image;
  typedef typename SceneGenerator::ImageContainer ImageContainer;
  typedef typename SceneGenerator::Point3D Point3D;
  typedef typename SceneGenerator::Point3DContainer Point3DContainer;

  typedef hs::sfm::KeysGenerator<Scalar, ImageDimension> KeysGenerator;
  typedef typename KeysGenerator::Keys Keys;
  typedef typename KeysGenerator::KeysContainer KeysContainer;

  typedef MultipleViewVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::IntrinsicParams MultipleViewIntrinsicParams;
  typedef typename VectorFunction::IntrinsicParamsContainer
                   MultipleViewIntrinsicParamsContainer;
  typedef typename VectorFunction::ExtrinsicParams MultipleViewExtrinsicParams;
  typedef typename VectorFunction::ExtrinsicParamsContainer
                   MultipleViewExtrinsicParamsContainer;

  MultipleViewSyntheicDataGenerator(
    Scalar focal_length_in_metre,
    size_t number_of_strps,
    size_t number_of_cameras_in_strip,
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    Scalar lateral_overlap_ratio,
    Scalar longitudinal_overlap_ratio,
    Scalar scene_max_height,
    Scalar camera_height_stddev,
    Scalar camera_plannar_stddev,
    Scalar camera_rot_stddev,
    Scalar north_west_angle)
    : scene_generator_(focal_length_in_metre,
                       number_of_strps,
                       number_of_cameras_in_strip,
                       ground_resolution,
                       image_width,
                       image_height,
                       pixel_size,
                       1,
                       lateral_overlap_ratio,
                       longitudinal_overlap_ratio,
                       scene_max_height,
                       camera_height_stddev,
                       camera_plannar_stddev,
                       camera_rot_stddev,
                       north_west_angle),
      keys_generator_(image_width, image_height){}

  Err operator()(VectorFunction& vector_function, XVector& x, YVector& y)
  {
    IntrinsicParamsContainer intrinsic_params_set;
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    Point3DContainer points;
    if (scene_generator_(intrinsic_params_set, extrinsic_params_set,
                         images, points) != 0)
    {
      return -1;
    }

    KeysContainer image_keys;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;
    if (keys_generator_(intrinsic_params_set, extrinsic_params_set, points,
                        image_keys, tracks, camera_views) != 0)
    {
      return -1;
    }

    MultipleViewIntrinsicParamsContainer multiple_view_intrinsic_params_set;
    MultipleViewExtrinsicParamsContainer multiple_view_extrinsic_params_set;
    for (size_t i = 0; i < intrinsic_params_set.size(); i++)
    {
      if (!camera_views[i].empty())
      {
        multiple_view_intrinsic_params_set.push_back(intrinsic_params_set[i]);
        multiple_view_extrinsic_params_set.push_back(extrinsic_params_set[i]);
      }
    }

    vector_function.set_intrinsic_params_set(
      multiple_view_intrinsic_params_set);
    vector_function.set_extrinsic_params_set(
      multiple_view_extrinsic_params_set);

    x = points[0];

    y.resize(Index(2 * multiple_view_intrinsic_params_set.size()));
    Index view_id = 0;
    for (size_t i = 0; i < intrinsic_params_set.size(); i++)
    {
      if (!camera_views[i].empty())
      {
        y.segment(view_id * 2, 2) = image_keys[i][0];
        view_id++;
      }
    }

    return 0;
  }

private:
  SceneGenerator scene_generator_;
  KeysGenerator keys_generator_;
};

}
}
}


#endif