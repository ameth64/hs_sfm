#ifndef _UNIT_TEST_TEST_MULTIPLE_VIEW_BASE_HPP_
#define _UNIT_TEST_TEST_MULTIPLE_VIEW_BASE_HPP_

#include <gtest/gtest.h>

#include "hs_sfm/synthetic/flight_generator.hpp"
#include "hs_sfm/synthetic/multiple_camera_keyset_generator.hpp"
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

  typedef hs::sfm::synthetic::FlightGenerator<Scalar, ImageDimension>
          FlightGenerator;
  typedef typename FlightGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename FlightGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename FlightGenerator::Image Image;
  typedef typename FlightGenerator::ImageContainer ImageContainer;
  typedef typename FlightGenerator::Point3D Point3D;
  typedef typename FlightGenerator::Point3DContainer Point3DContainer;

  typedef hs::sfm::synthetic::MultipleCameraKeysetGenerator<Scalar,
                                                            ImageDimension>
          KeysetGenerator;
  typedef typename KeysetGenerator::IntrinsicParams IntrinsicParams;
  typedef typename KeysetGenerator::IntrinsicParamsContainer
          IntrinsicParamsContainer;
  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

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
    Scalar north_west_angle,
    const IntrinsicParamsContainer& intrinsic_params_set)
    : flight_generator_(focal_length_in_metre,
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
                        camera_rot_stddev)
    , north_west_angle_(north_west_angle)
    , intrinsic_params_set_(intrinsic_params_set){}

  Err operator()(VectorFunction& vector_function, XVector& x, YVector& y)
  {
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    flight_generator_.GenerateCameras(extrinsic_params_set, images);
    while (1)
    {
      Point3DContainer points;
      flight_generator_.GeneratePoints(1, points);
      std::vector<size_t> image_intrinsic_map;
      for (size_t i = 0; i < extrinsic_params_set.size(); i++)
      {
        image_intrinsic_map.push_back(i);
      }
      KeysetContainer keysets;
      TrackContainer tracks;
      CameraViewContainer camera_views;
      if (keys_generator_(intrinsic_params_set_,
                          extrinsic_params_set,
                          images,
                          points,
                          image_intrinsic_map,
                          keysets,
                          tracks,
                          camera_views) != 0) continue;
      if (tracks[0].size() == extrinsic_params_set.size())
      {
        vector_function.set_intrinsic_params_set(intrinsic_params_set_);
        vector_function.set_extrinsic_params_set(extrinsic_params_set);
        y.resize(Index(extrinsic_params_set.size() * 2));
        for (size_t i = 0; i < extrinsic_params_set.size(); i++)
        {
          y[Index(i * 2) + 0] = keysets[i][0][0];
          y[Index(i * 2) + 1] = keysets[i][0][1];
        }
        x = points[0];
        break;
      }
    }
    return 0;
  }

private:
  FlightGenerator flight_generator_;
  KeysetGenerator keys_generator_;
  Scalar north_west_angle_;
  IntrinsicParamsContainer intrinsic_params_set_;
};

}
}
}


#endif
