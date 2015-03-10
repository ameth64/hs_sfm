#ifndef _HS_SFM_UNIT_TEST_ESSENTIAL_SYNTHETIC_DATA_GENERATOR_HPP_
#define _HS_SFM_UNIT_TEST_ESSENTIAL_SYNTHETIC_DATA_GENERATOR_HPP_

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/synthetic/flight_generator.hpp"
#include "hs_sfm/synthetic/multiple_camera_keyset_generator.hpp"

namespace hs
{
namespace sfm
{
namespace essential
{

template <typename _Scalar, typename _ImageDimension>
class SyntheticDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;

private:
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
  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;
public:
  typedef typename KeysetGenerator::IntrinsicParams IntrinsicParams;
  typedef typename KeysetGenerator::IntrinsicParamsContainer
    IntrinsicParamsContainer;
  typedef typename IntrinsicParams::KMatrix KMatrix;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef std::pair<Key, Key> KeyPair;
  typedef EIGEN_STD_VECTOR(KeyPair) KeyPairContainer;

  SyntheticDataGenerator(
    Scalar focal_length_in_metre,
    size_t number_of_strps,
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    size_t number_of_points,
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
                        number_of_strps > 1 ? 1 : 2,
                        ground_resolution,
                        image_width,
                        image_height,
                        pixel_size,
                        number_of_points,
                        lateral_overlap_ratio,
                        longitudinal_overlap_ratio,
                        scene_max_height,
                        camera_height_stddev,
                        camera_plannar_stddev,
                        camera_rot_stddev),
      intrinsic_params_set_(intrinsic_params_set) {}

  int operator() (KeyPairContainer& key_pairs,
                  KMatrix& K_left, KMatrix& K_right)
  {
    //生成相机参数和三维点
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    Point3DContainer points;
    if (flight_generator_(extrinsic_params_set, images, points) != 0)
    {
      return -1;
    }

    //生成特征点
    std::vector<size_t> image_intrinsic_map;
    image_intrinsic_map.push_back(0);
    image_intrinsic_map.push_back(1);
    KeysetContainer keysets;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;
    if (keys_generator_(intrinsic_params_set_,
                        extrinsic_params_set,
                        images,
                        points,
                        image_intrinsic_map,
                        keysets,
                        tracks,
                        camera_views) != 0)
    {
      std::cout<<"keys generator failed!\n";
      return -1;
    }

    EIGEN_MATRIX(Scalar, 2, 2) key_covariance;
    key_covariance.setIdentity();
    for (size_t i = 0; i < tracks.size(); i++)
    {
      if (tracks[i].size() == 2)
      {
        KeyPair key_pair;
        size_t key_left_id = tracks[i][0].second;
        size_t key_right_id = tracks[i][1].second;
        EIGEN_VECTOR(Scalar, 2) noised_key;
        hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
          keysets[0][key_left_id], key_covariance, noised_key);
        key_pair.first = noised_key;
        hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
          keysets[1][key_right_id], key_covariance, noised_key);
        key_pair.second = noised_key;
        key_pairs.push_back(key_pair);
      }
    }

    K_left = intrinsic_params_set_[0].GetKMatrix();
    K_right = intrinsic_params_set_[1].GetKMatrix();

    return 0;
  }

private:
  FlightGenerator flight_generator_;
  KeysetGenerator keys_generator_;
  IntrinsicParamsContainer intrinsic_params_set_;
};

}
}
}

#endif