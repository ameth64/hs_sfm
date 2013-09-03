#ifndef _UNIT_TEST_TEST_MULTIPLE_VIEW_BASE_HPP_
#define _UNIT_TEST_TEST_MULTIPLE_VIEW_BASE_HPP_

#include <gtest/gtest.h>

#include "hs_sfm/utility/synthetic_scene_generator.hpp"
#include "hs_sfm/triangulate/multiple_view_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace triangulate
{

template <typename _Scalar, typename _ImgDim>
class MultipleViewSyntheicDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImgDim ImgDim;

  typedef int Err;

  typedef hs::sfm::SceneGenerator<Scalar, ImgDim> SceneGenerator;
  typedef typename SceneGenerator::Intrin Intrin;
  typedef typename SceneGenerator::IntrinContainer IntrinContainer;
  typedef typename SceneGenerator::Extrin Extrin;
  typedef typename SceneGenerator::ExtrinContainer ExtrinContainer;
  typedef typename SceneGenerator::Image Image;
  typedef typename SceneGenerator::ImageContainer ImageContainer;
  typedef typename SceneGenerator::Pt3D Pt3D;
  typedef typename SceneGenerator::Pt3DContainer Pt3DContainer;

  typedef hs::sfm::KeysGenerator<Scalar, ImgDim> KeysGenerator;
  typedef typename KeysGenerator::ImgKeys ImgKeys;
  typedef typename KeysGenerator::ImgKeysContainer ImgKeysContainer;
  typedef typename KeysGenerator::Track Track;
  typedef typename KeysGenerator::TrackContainer TrackContainer;
  typedef typename KeysGenerator::CamViewContainer CamViewContainer;

  typedef MultipleViewVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVec XVector;
  typedef typename VectorFunction::YVec YVector;
  typedef typename VectorFunction::CameraIntrin MultipleViewIntrin;
  typedef typename VectorFunction::CameraIntrinContainer
                   MultipleViewIntrinContainer;
  typedef typename VectorFunction::CameraExtrin MultipleViewExtrin;
  typedef typename VectorFunction::CameraExtrinContainer
                   MultipleViewExtrinContainer;

  MultipleViewSyntheicDataGenerator(
    Scalar f, size_t strip_num,
    size_t cams_num_in_strip,
    Scalar ground_resolution,
    ImgDim img_width, ImgDim img_height,
    Scalar pixel_size,
    Scalar lateral_overlap,
    Scalar longitudinal_overlap,
    Scalar scene_max_height,
    Scalar cam_height_stddev,
    Scalar cam_plannar_stddev,
    Scalar cam_rot_stddev,
    Scalar nwAngle)
    : scene_generator_(f, strip_num, cams_num_in_strip, ground_resolution,
                       img_width, img_height, pixel_size, 1,
                       lateral_overlap, longitudinal_overlap,
                       scene_max_height, cam_height_stddev, cam_plannar_stddev,
                       cam_rot_stddev, nwAngle),
      keys_generator_(img_width, img_height){}

  Err operator()(VectorFunction& vector_function, XVector& x, YVector& y)
  {
    IntrinContainer intrins;
    ExtrinContainer extrins;
    ImageContainer images;
    Pt3DContainer pts;
    if (scene_generator_(intrins, extrins, images, pts) != 0)
    {
      return -1;
    }

    ImgKeysContainer img_keys;
    TrackContainer tracks;
    CamViewContainer cam_views;
    if (keys_generator_(intrins, extrins, pts,
                        img_keys, tracks, cam_views) != 0)
    {
      return -1;
    }

    MultipleViewIntrinContainer multiple_view_intrins;
    MultipleViewExtrinContainer multiple_view_extrins;
    for (size_t i = 0; i < intrins.size(); i++)
    {
      if (!cam_views[i].empty())
      {
        multiple_view_intrins.push_back(intrins[i]);
        multiple_view_extrins.push_back(extrins[i]);
      }
    }

    vector_function.SetIntrins(multiple_view_intrins);
    vector_function.SetExtrins(multiple_view_extrins);

    x = pts[0];

    y.resize(Index(2 * multiple_view_intrins.size()));
    Index view_id = 0;
    for (size_t i = 0; i < intrins.size(); i++)
    {
      if (!cam_views[i].empty())
      {
        y.segment(view_id * 2, 2) = img_keys[i][0];
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