#ifndef _HS_SFM_SYNTHETIC_MULTIPLE_CAMERA_KEYSET_GENERATOR_HPP_
#define _HS_SFM_SYNTHETIC_MULTIPLE_CAMERA_KEYSET_GENERATOR_HPP_

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/radial_distortor.hpp"
#include "hs_sfm/sfm_utility/decentering_distortor.hpp"
#include "hs_sfm/sfm_utility/projective_functions.hpp"

namespace hs
{
namespace sfm
{
namespace synthetic
{

template <typename _Scalar, typename _ImageDimension>
class MultipleCameraKeysetGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef ImageParams<ImageDimension> Image;
  typedef EIGEN_STD_VECTOR(Image) ImageContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_STD_VECTOR(Point3D) Point3DContianer;
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;

private:
  typedef typename Keyset::Key Key;
  typedef typename Keyset::KeysContainer KeyContainer;
  typedef ProjectiveFunctions<Scalar> ProjectiveFunctions;

public:
  Err operator() (const IntrinsicParamsContainer& intrinsic_params_set,
                  const ExtrinsicParamsContainer& extrinsic_params_set,
                  const ImageContainer& images,
                  const Point3DContianer& points,
                  const std::vector<size_t>& image_intrinsic_map,
                  KeysetContainer& keysets,
                  TrackContainer& tracks,
                  CameraViewContainer& camera_views) const
  {
    size_t number_of_images = images.size();
    if (number_of_images != extrinsic_params_set.size() ||
        number_of_images != image_intrinsic_map.size())
    {
      return -1;
    }

    EIGEN_STD_VECTOR(KeyContainer) keys_vector;
    keys_vector.resize(number_of_images);	//每张图像使用一个KeyContainer保存关键点
    size_t number_of_points = points.size();
    tracks.resize(number_of_points);
    camera_views.resize(number_of_images);
    for (size_t i = 0; i < number_of_points; i++)
    {
      const Point3D& point = points[i];
      for (size_t j = 0; j < number_of_images; j++)
      {
        size_t intrinsic_id = image_intrinsic_map[j];
        const IntrinsicParams& intrinsic_params =
          intrinsic_params_set[intrinsic_id];
        const ExtrinsicParams& extrinsic_params =
          extrinsic_params_set[j];
        const Image& image = images[j];
        Key key =
          ProjectiveFunctions::WorldPointProjectToImageKeyNoDistort(
            intrinsic_params,
            extrinsic_params,
            point);	//以不加扭曲的方式将空间点投影至像平面并判断是否在成像区域内

        if (key[0] > 0 && key[0] < Scalar(image.m_width) &&
            key[1] > 0 && key[1] < Scalar(image.m_height))
        {
          key = ProjectiveFunctions::WorldPointProjectToImageKey(
                  intrinsic_params,
                  extrinsic_params,
                  point);
          if (key[0] > 0 && key[0] < Scalar(image.m_width) &&
              key[1] > 0 && key[1] < Scalar(image.m_height))
          {
            tracks[i].push_back(
              std::make_pair(j, keys_vector[j].size()));	//记录图像编号 - key点在该图像的特征点列表中的序号, 这两个变量pair也就是track
            camera_views[j].push_back(
              std::make_pair(i, keys_vector[j].size()));	//记录3D点编号 - key点序号
            keys_vector[j].push_back(key);	//第j张图像的KeyContainer记录key点
          }
        }
      }// for (size_t j = 0; j < number_of_images; j++)
    }// for (size_t i = 0; i < number_of_points; i++)

    keysets.clear();
    for (size_t i = 0; i < number_of_images; i++)
    {
      keysets.push_back(Keyset(keys_vector[i]));
    }

    return 0;
  }
};

}
}
}

#endif
