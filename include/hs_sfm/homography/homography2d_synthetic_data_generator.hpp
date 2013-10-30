#ifndef _HS_SFM_HOMOGRAPHY_HOMOGRAPHY2D_SYNTHETIC_DATA_GENERATOR_HPP_
#define _HS_SFM_HOMOGRAPHY_HOMOGRAPHY2D_SYNTHETIC_DATA_GENERATOR_HPP_

#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/synthetic/keyset_generator.hpp"
#include "hs_sfm/homography/homography2d_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace homography
{

template <typename _Scalar, typename _ImageDimension>
class Homography2DSyntheticDataGenerator
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

  typedef Homography2DVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::HomogeneousMatrix HomogeneousMatrix;
  typedef typename VectorFunction::Key Key;
  typedef typename VectorFunction::HomogeneousKey HomogeneousKey;

public:
  Homography2DSyntheticDataGenerator(
    Scalar focal_length_in_metre,
    size_t number_of_strips,
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    Scalar lateral_overlap_ratio,
    Scalar longitudinal_overlap_ratio,
    Scalar scene_max_height,
    Scalar camera_height_stddev,
    Scalar camera_planar_stddev,
    Scalar camera_rotation_stddev,
    Scalar north_west_angle,
    size_t number_of_keys,
    Scalar points_height)
    : scene_generator_(focal_length_in_metre,
                       number_of_strips,
                       number_of_strips > 1 ? 1 : 2,
                       ground_resolution,
                       image_width,
                       image_height,
                       pixel_size,
                       1,
                       lateral_overlap_ratio,
                       longitudinal_overlap_ratio,
                       scene_max_height,
                       camera_height_stddev,
                       camera_planar_stddev,
                       camera_rotation_stddev,
                       north_west_angle),
      keys_generator_(image_width, image_height),
      number_of_keys_(number_of_keys),
      points_height_(points_height) {}

public:
  Err operator()(VectorFunction& vector_function,
                 XVector& x,
                 YVector& y) const
  {
    IntrinsicParamsContainer intrinsic_params_set;
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    scene_generator_.GenerateSceneCameras(intrinsic_params_set,
                                          extrinsic_params_set,
                                          images);
    Point3DContainer points;
    scene_generator_.GenerateScenePlanarPoints(number_of_keys_,
                                               points_height_,
                                               points);

    KeysetContainer keysets;
    TrackContainer tracks;
    CameraViewContainer camera_views;
    if (keys_generator_(intrinsic_params_set,
                        extrinsic_params_set,
                        points,
                        keysets,
                        tracks,
                        camera_views) != 0)
    {
      return -1;
    }

    typename TrackContainer::const_iterator itr_track = tracks.begin();
    typename TrackContainer::const_iterator itr_track_end = tracks.end();
    std::vector<size_t> point_maps(points.size(), 0);
    size_t number_of_valid_track = 0;
    for (size_t i = 0; itr_track != itr_track_end; ++itr_track)
    {
      if (itr_track->size() == 2)
      {
        point_maps[i] = number_of_valid_track;
        number_of_valid_track++;
        i++;
      }
    }

    vector_function.set_number_of_keys(number_of_valid_track);
    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    x.resize(x_size);
    y.resize(y_size);

    typename IntrinsicParams::KMatrix K0 =
      intrinsic_params_set[0].GetKMatrix();
    typename IntrinsicParams::KMatrix K1 =
      intrinsic_params_set[1].GetKMatrix();
    EIGEN_MATRIX(Scalar, 3, 3) R0 =
      extrinsic_params_set[0].rotation();
    EIGEN_MATRIX(Scalar, 3, 3) R1 =
      extrinsic_params_set[1].rotation();
    EIGEN_VECTOR(Scalar, 3) c0 =
      extrinsic_params_set[0].position();
    EIGEN_VECTOR(Scalar, 3) c1 =
      extrinsic_params_set[1].position();
    EIGEN_VECTOR(Scalar, 3) t0 = -R0 * c0;
    EIGEN_VECTOR(Scalar, 3) t1 = -R1 * c1;
    R0.col(2) += t0 / points_height_;
    R1.col(2) += t1 / points_height_;

    HomogeneousMatrix H0 = K0 * R0;
    HomogeneousMatrix H1 = K1 * R1;
    HomogeneousMatrix H = H1 * H0.inverse();
    vector_function.SetHomographyMatrix(H, x);

    itr_track = tracks.begin();
    itr_track_end = tracks.end();
    for (Index i = 0; itr_track != itr_track_end; ++itr_track)
    {
      if (itr_track->size() == 2)
      {
        x.template segment<2>(i * 2) = keysets[0][(*itr_track)[0].second];
        y.template segment<2>(i * 4) = keysets[0][(*itr_track)[0].second];
        y.template segment<2>(i * 4 + 2) = keysets[1][(*itr_track)[1].second];

        i++;
      }
    }
 
    return 0;
  }

private:
  SceneGenerator scene_generator_;
  KeysetGenerator keys_generator_;
  size_t number_of_keys_;
  Scalar points_height_;
};

}
}
}

#endif