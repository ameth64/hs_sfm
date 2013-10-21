#include <gtest/gtest.h>

#include "hs_sfm/sfm_utility/synthetic_scene_generator.hpp"
#include "hs_sfm/sfm_utility/sfm_file_io.hpp"

namespace
{

template <typename _Scalar>
struct PointsProjectToCameraResidual
{
  typedef _Scalar Scalar;
  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef hs::sfm::ImageKeys<Scalar> ImgKeys;
  typedef EIGEN_STD_VECTOR(ImgKeys) ImgKeysContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef Point3D KeyHomogeneous;
  typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;
  typedef std::vector<std::pair<size_t, size_t> > Track;
  typedef std::vector<Track> TrackContainer;
  typedef hs::sfm::CameraFunctions<Scalar> Camera;
  typedef typename Camera::ProjectionMatrix ProjectionMatrix;

  Scalar operator()(const IntrinsicParamsContainer& intrinsic_params_set,
                    const ExtrinsicParamsContainer& extrinsic_params_set,
                    const Point3DContainer& points,
                    const ImgKeysContainer& img_keys_set,
                    const TrackContainer& tracks) const
  {
    size_t number_of_points = points.size();
    size_t number_of_cameras = intrinsic_params_set.size();
    if (intrinsic_params_set.size() != number_of_cameras || img_keys_set.size() != number_of_cameras ||
      tracks.size() != number_of_points)
    {
      return Scalar(-1);
    }

    Scalar distance_sum = Scalar(0);
    size_t number_of_measures = 0;
    for (size_t i = 0; i < number_of_points; i++)
    {
      size_t number_of_views = tracks[i].size();
      number_of_measures += number_of_views;
      for (size_t j = 0; j < number_of_views; j++)
      {
        size_t img_id = tracks[i][j].first;
        size_t key_id = tracks[i][j].second;

        Key key = img_keys_set[img_id][key_id];

        ProjectionMatrix P =
          Camera::GetProjectionMatrix(intrinsic_params_set[img_id],
                                      extrinsic_params_set[img_id]);
        KeyHomogeneous key_homogeneous = P.block(0, 0, 3, 3) * points[i] + 
                                         P.block(0, 3, 3, 1);
        key_homogeneous /= key_homogeneous(2);

        Scalar distance = (key_homogeneous.segment(0, 2) - key).squaredNorm();
        distance_sum += distance;
      }
    }

    Scalar residual =
      sqrt(distance_sum / Scalar(number_of_measures) / Scalar(2));

    return residual;
  }
};

TEST(TestSyntheticSceneGenerator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef hs::sfm::SceneGenerator<Scalar, ImageDimension> SceneGen;
  typedef hs::sfm::KeysGenerator<Scalar, ImageDimension> KeysGen;
  typedef hs::sfm::MatchGenerator MatchesGen;
  typedef hs::sfm::NoiseImageKeys<Scalar> NoiseImgKeys;
  typedef hs::sfm::NoiseExtrin<Scalar> NoiseExtr;

  typedef SceneGen::IntrinsicParams IntrinsicParams;
  typedef SceneGen::ExtrinsicParams ExtrinsicParams;
  typedef SceneGen::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef SceneGen::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef SceneGen::ImageContainer ImageContainer;
  typedef SceneGen::Point3DContainer Point3DContainer;

  typedef KeysGen::Keys ImgKeys;
  typedef KeysGen::KeysContainer ImgKeysContainer;
  typedef KeysGen::TrackContainer TrackContainer;
  typedef KeysGen::CameraViewContainer CameraViewContainer;

  typedef hs::sfm::MatchContainer MatchContainer;

  typedef NoiseImgKeys::Covariance KeyCovariance;

  typedef NoiseExtr::Covariance ExtrinCovariance;

  typedef std::pair<size_t, size_t> CameraPairKey;
  typedef EIGEN_STD_MAP(CameraPairKey, ExtrinsicParams) CameraPairContainer;

  typedef PointsProjectToCameraResidual<Scalar> Residual;

  Scalar focal_length_in_metre = 0.006;
  size_t number_of_strips = 15;
  size_t number_of_cameras_in_strip = 20;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 3648;
  ImageDimension image_height = 2736;
  Scalar pixel_size = 0.00000203311408298266;
  size_t number_of_points = 2000;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 2;
  Scalar camera_planar_stddev = 2;
  Scalar camera_rotation_stddev = 10;
  Scalar north_west_angle = 60;

  SceneGen scene_generator(focal_length_in_metre,
                           number_of_strips,
                           number_of_cameras_in_strip,
                           ground_resolution, 
                           image_width,
                           image_height,
                           pixel_size,
                           number_of_points,
                           lateral_overlap_ratio,
                           longitudinal_overlap_ratio,
                           scene_max_height,
                           camera_height_stddev,
                           camera_planar_stddev,
                           camera_rotation_stddev,
                           north_west_angle);

  IntrinsicParamsContainer intrinsic_params_set;
  ExtrinsicParamsContainer extrinsic_params_set;
  ImageContainer images;
  Point3DContainer points;
  scene_generator(intrinsic_params_set, extrinsic_params_set, images, points);
  //保存文件方便可视化观察
  //sfm::SaveXugFile<Scalar, ImageDimension> saver;
  //saver("TestSyntheticSceneGenerator/SimpleTest/test.xug",
  //    intrinsic_params_set, extrinsic_params_set, images, points, 10);

  KeysGen keys_generator(image_width, image_height);
  ImgKeysContainer image_keys_set;
  TrackContainer tracks;
  CameraViewContainer camera_views;
  keys_generator(intrinsic_params_set,
                 extrinsic_params_set, points,
                 image_keys_set,
                 tracks,
                 camera_views);

  MatchesGen matches_generator;
  MatchContainer matches;
  matches_generator(tracks, matches);

  //特征点加入噪声
  ImgKeysContainer noised_image_keys_set = image_keys_set;
  size_t number_of_cameras = noised_image_keys_set.size();
  KeyCovariance key_covariance = KeyCovariance::Identity();
  NoiseImgKeys noise_image_keys(key_covariance);
  for (size_t i = 0; i < number_of_cameras; i++)
  {
    noise_image_keys(noised_image_keys_set[i]);
  }

  //外参数加入噪声，作为POS信息
  ExtrinsicParamsContainer prior_extrinsic_params_set = extrinsic_params_set;
  ExtrinCovariance rotation_covariance = ExtrinCovariance::Identity();
  rotation_covariance *= Scalar(5) / 180 * Scalar(M_PI);
  ExtrinCovariance position_covariance = ExtrinCovariance::Identity();
  position_covariance *= Scalar(10);
  NoiseExtr noise_extr(rotation_covariance, position_covariance);
  for (size_t i = 0; i < number_of_cameras; i++)
  {
    noise_extr(prior_extrinsic_params_set[i]);
  }

  //保存文件方便可视化观察
  //saver("TestSyntheticSceneGenerator/SimpleTest/test_noise.xug",
  //  intrinsic_params_set, prior_extrinsic_params_set, images, points, 10);

}

}
