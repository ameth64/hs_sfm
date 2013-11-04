#ifndef _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTAINED_SYNTHETIC_GENERATOR_HPP_
#define _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTAINED_SYNTHETIC_GENERATOR_HPP_

#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/synthetic/keyset_generator.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

template <typename _Scalar, typename _ImageDimension>
class IntrinsicConstrainedSyntheticGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

  typedef hs::sfm::synthetic::SceneGenerator<Scalar, ImageDimension>
          SceneGenerator;
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

  typedef hs::sfm::synthetic::KeysetGenerator<Scalar, ImageDimension>
          KeysetGenerator;
  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

  typedef IntrinsicConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;

public:
  IntrinsicConstrainedSyntheticGenerator(
    Scalar focal_length_in_metre,
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    size_t number_of_points,
    Scalar scene_max_height,
    Scalar camera_height_stddev,
    Scalar camera_planar_stddev,
    Scalar camera_rot_stddev,
    Scalar north_west_angle,
    Scalar skew,
    Scalar principal_point_x,
    Scalar principal_point_y,
    Scalar pixel_ratio)
    : scene_generator_(focal_length_in_metre,
                       1, 1,
                       ground_resolution,
                       image_width,
                       image_height,
                       pixel_size,
                       number_of_points,
                       0.6, 0.8,
                       scene_max_height,
                       camera_height_stddev,
                       camera_planar_stddev,
                       camera_rot_stddev,
                       north_west_angle),
      keyset_generator_(image_width, image_height),
      skew_(skew),
      principal_point_x_(principal_point_x),
      principal_point_y_(principal_point_y),
      pixel_ratio_(pixel_ratio){}

  Err operator() (VectorFunction& vector_function,
                  XVector& x,
                  YVector& y) const
  {
    IntrinsicParamsContainer intrinsic_params_set;
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    Point3DContainer points;
    if (scene_generator_(intrinsic_params_set,
                         extrinsic_params_set,
                         images,
                         points) != 0) return -1;

    intrinsic_params_set[0].set_skew(skew_);
    intrinsic_params_set[0].set_principal_point_x(principal_point_x_);
    intrinsic_params_set[0].set_principal_point_y(principal_point_y_);
    intrinsic_params_set[0].set_pixel_ratio(pixel_ratio_);

    KeysetContainer keysets;
    TrackContainer tracks;
    CameraViewContainer camera_views;
    if (keyset_generator_(intrinsic_params_set,
                          extrinsic_params_set,
                          points, keysets, tracks, camera_views) != 0)
      return -1;

    vector_function.points() = points;
    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    x.resize(x_size);
    y.resize(y_size);
    x[0] = extrinsic_params_set[0].rotation()[0];
    x[1] = extrinsic_params_set[0].rotation()[1];
    x[2] = extrinsic_params_set[0].rotation()[2];
    x[3] = extrinsic_params_set[0].position()[0];
    x[4] = extrinsic_params_set[0].position()[1];
    x[5] = extrinsic_params_set[0].position()[2];
    x[6] = intrinsic_params_set[0].focal_length();
    x[7] = intrinsic_params_set[0].skew();
    x[8] = intrinsic_params_set[0].principal_point_x();
    x[9] = intrinsic_params_set[0].principal_point_y();
    x[10] = intrinsic_params_set[0].pixel_ratio();

    Index number_of_keys = Index(camera_views[0].size());
    y[number_of_keys * 2 + 0] = intrinsic_params_set[0].focal_length();
    y[number_of_keys * 2 + 1] = intrinsic_params_set[0].skew();
    y[number_of_keys * 2 + 2] = intrinsic_params_set[0].principal_point_x();
    y[number_of_keys * 2 + 3] = intrinsic_params_set[0].principal_point_y();
    y[number_of_keys * 2 + 4] = intrinsic_params_set[0].pixel_ratio();
    for (Index i = 0; i < number_of_keys; i++)
    {
      y[i * 2 + 0] = keysets[0][i][0];
      y[i * 2 + 1] = keysets[0][i][1];
    }

    return 0;
  }

private:
  SceneGenerator scene_generator_;
  KeysetGenerator keyset_generator_;
  Scalar skew_;
  Scalar principal_point_x_;
  Scalar principal_point_y_;
  Scalar pixel_ratio_;
};

}
}
}

#endif
