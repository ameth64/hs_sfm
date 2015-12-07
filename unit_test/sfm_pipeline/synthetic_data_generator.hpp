#ifndef _HS_SFM_UNIT_TEST_SFM_PIPELINE_SYNTHETIC_DATA_GENERATOR_HPP_
#define _HS_SFM_UNIT_TEST_SFM_PIPELINE_SYNTHETIC_DATA_GENERATOR_HPP_

#include <string>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"
#include "hs_math/geometry/euler_angles.hpp"

#include "hs_sfm/synthetic/multiple_camera_keyset_generator.hpp"
//#include "hs_sfm/synthetic/relative_generator.hpp"
#include "hs_sfm/synthetic/multiple_flight_generator.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"
#if 1
#include "hs_sfm/sfm_file_io/scene_ply_saver.hpp"
#endif

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
struct RichTrack	//在Track的基础上封装的数据结构
{
  typedef _Scalar Scalar;
  hs::sfm::Track track;
  EIGEN_VECTOR(Scalar, 3) point;

  bool operator < (const RichTrack<Scalar>& other) const
  {
    return (track < other.track);
  }
};

/************************************************************************/
/* 
航线数据生成类, 下辖MultipleFlightGenerator, FlightGenerator等子类

*/
/************************************************************************/
template <typename _Scalar, typename _ImageDimension>
class SyntheticDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
  typedef hs::sfm::synthetic::MultipleFlightGenerator<Scalar, ImageDimension>
          SceneGenerator;
  typedef hs::sfm::synthetic::MultipleCameraKeysetGenerator<Scalar,
                                                            ImageDimension>
          KeysetGenerator;
  typedef RichTrack<Scalar> RichTrackType;
  typedef EIGEN_STD_VECTOR(RichTrackType) RichTrackContainer;
  //typedef hs::sfm::synthetic::RelativeGenerator<Scalar> RelativeGenerator;

public:
  typedef typename KeysetGenerator::IntrinsicParams IntrinsicParams;
  typedef typename IntrinsicParams::KMatrix KMatrix;
  typedef typename KeysetGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename SceneGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename SceneGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename SceneGenerator::Image Image;
  typedef typename SceneGenerator::ImageContainer ImageContainer;
  typedef typename SceneGenerator::Point3D Point3D;
  typedef typename SceneGenerator::Point3DContainer Point3DContainer;
  typedef typename SceneGenerator::FlightGenerator FlightGenerator;
  typedef typename SceneGenerator::FlightGeneratorContainer
                   FlightGeneratorContainer;

  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

  //typedef typename RelativeGenerator::RMatrix RMatrix;
  //typedef typename RelativeGenerator::Translate Translate;

public:
  SyntheticDataGenerator(
    Scalar flight_longitudinal_overlap_ratio,
    Scalar flight_lateral_overlap_ratio,
    Scalar north_west_angle,
    Scalar north_west_angle_stddev,
    Scalar offset_stddev,
    const FlightGeneratorContainer& flight_generators,
    Scalar outlier_ratio,
    Scalar key_stddev,
    const IntrinsicParamsContainer& intrinsic_params_set_true,
    size_t number_of_points,
    size_t number_of_gcps,
    size_t number_of_check_points)
  : scene_generator_(flight_longitudinal_overlap_ratio,
                     flight_lateral_overlap_ratio,
                     north_west_angle,
                     north_west_angle_stddev,
                     offset_stddev,
                     flight_generators),
    outlier_ratio_(outlier_ratio),
    key_stddev_(key_stddev),
    intrinsic_params_set_true_(intrinsic_params_set_true),
    number_of_points_(number_of_points),
    number_of_gcps_(number_of_gcps),
    number_of_check_points_(number_of_check_points) {}

public:
  Err Generate(std::vector<size_t>& image_intrinsic_map,
               hs::sfm::MatchContainer& matches,
               KeysetContainer& keysets_noised,
               Point3DContainer& gcps,
               hs::sfm::TrackContainer& tracks_gcp,
               KeysetContainer& keysets_gcp_noised,
               Point3DContainer& check_points,
               hs::sfm::TrackContainer& tracks_check_point,
               KeysetContainer& keysets_checkout_point_noised,
               ExtrinsicParamsContainer& extrinsic_params_set_absolute_true,
               Point3DContainer& points_absolute_true) const
  {
    ImageContainer images;
    Point3DContainer points_absolute;
    KeysetContainer keysets_true;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;

    Err result = 0;
    while (1)
    {
      result = GenerateAbsoluteScene(extrinsic_params_set_absolute_true,
                                     images,
                                     points_absolute,
                                     keysets_true,
                                     image_intrinsic_map,
                                     tracks,
                                     camera_views);	//生成空间3D点及其在每个图像上的同名点等基础数据
      if (result != 0) break;


      RichTrackContainer rich_tracks;
      for (size_t i = 0; i < number_of_points_; i++)	//构造RichTrack集合
      {
        if (tracks[i].size() > 1)
        {
          RichTrackType rich_track;
          rich_track.point = points_absolute[i];
          rich_track.track = tracks[i];
          std::sort(rich_track.track.begin(), rich_track.track.end());	//对track的成员排序
          rich_tracks.push_back(rich_track);
        }
      }
      std::sort(rich_tracks.begin(), rich_tracks.end());	
      size_t number_of_rich_tracks = rich_tracks.size();
      TrackContainer tracks_true(number_of_rich_tracks);

      points_absolute_true.resize(number_of_rich_tracks);
      for (size_t i = 0; i < number_of_rich_tracks; i++)
      {
        tracks_true[i] = rich_tracks[i].track;
        points_absolute_true[i] = rich_tracks[i].point;
      }

      hs::sfm::MatchesTracksConvertor matches_tracks_convertor;
      result = matches_tracks_convertor(tracks_true, matches);	//从TrackContainer构造MatchesTracksConvertor,
      if (result != 0) break;

      result = GenerateNoisedKeysets(keysets_true, images, keysets_noised);	//在生成的真值点的基础上添加随机噪声用于模拟
      if (result != 0) break;

      KeysetContainer keysets_gcp_true;
      result = GenerateGCPData(intrinsic_params_set_true_,
                               extrinsic_params_set_absolute_true,
                               images,
                               image_intrinsic_map,
                               number_of_gcps_,
                               keysets_gcp_true,
                               gcps,
                               tracks_gcp);	//生成地面像控点数据
      if (result != 0) break;

      result = GenerateNoisedKeysets(keysets_gcp_true, images,
                                     keysets_gcp_noised, false);
      if (result != 0) break;

      KeysetContainer keysets_check_point_true;
      result = GenerateGCPData(intrinsic_params_set_true_,
                               extrinsic_params_set_absolute_true,
                               images,
                               image_intrinsic_map,
                               number_of_check_points_,
                               keysets_check_point_true,
                               check_points,
                               tracks_check_point);	
      if (result != 0) break;

      result = GenerateNoisedKeysets(keysets_check_point_true, images,
                                     keysets_checkout_point_noised, false);
      if (result != 0) break;

      break;
    }

    return result;
  }

private:
  Err GenerateAbsoluteScene(
    ExtrinsicParamsContainer& extrinsic_params_set_absolute,
    ImageContainer& images,
    Point3DContainer& points_absolute,
    KeysetContainer& keysets,
    std::vector<size_t>& image_intrinsic_map,
    hs::sfm::TrackContainer& tracks,
    hs::sfm::CameraViewContainer& camera_views) const	//该方法完成空间随机点的生成, 航线上各图像的相机参数矩阵计算以及空间点在这些图像上的投影, 并生成track等数据结构
  {
    if (scene_generator_.GeneratePoints(number_of_points_,
                                        points_absolute) != 0)	//生成空间随机点,注意并未使用FlightGenerator对象的方法
    {
      return -1;
    }

    size_t number_of_flights = scene_generator_.GetNumberOfFlights();
    image_intrinsic_map.clear();
    extrinsic_params_set_absolute.clear();
#if 1
    IntrinsicParamsContainer intrinsic_params_set_display;
#endif
    for (size_t i = 0; i < number_of_flights; i++)
    {
      ExtrinsicParamsContainer extrinsic_params_set_flight;
      ImageContainer images_flight;
      if (scene_generator_.GenerateExtrinsicParamsContainer(
            i, extrinsic_params_set_flight, images_flight))	//计算相机参数并映射至世界坐标系
      {
        return -1;
      }
      for (size_t j = 0; j < extrinsic_params_set_flight.size(); j++)
      {
        extrinsic_params_set_absolute.push_back(extrinsic_params_set_flight[j]);
        images.push_back(images_flight[j]);
        image_intrinsic_map.push_back(i);
#if 1
        intrinsic_params_set_display.push_back(intrinsic_params_set_true_[i]);
#endif
      }
    }

#if 1
    typedef hs::sfm::fileio::ScenePLYSaver<Scalar, ImageDimension> Saver;
    Saver saver(30);
    saver("incremental_synthetic.ply",
          intrinsic_params_set_display, extrinsic_params_set_absolute,
          images, points_absolute);
#endif

    if (keyset_generator_(intrinsic_params_set_true_,
                          extrinsic_params_set_absolute,
                          images,
                          points_absolute,
                          image_intrinsic_map,
                          keysets,
                          tracks,
                          camera_views) != 0) //将空间点映射至每个图像的像平面, 并生成对应的track等数据结构对象
    {
      return -1;
    }

    return 0;
  }

  Err GenerateGCPData(
    const IntrinsicParamsContainer& intrinsic_params_set,
    const ExtrinsicParamsContainer& extrinsic_params_set_absolute,
    const ImageContainer& images,
    const std::vector<size_t>& image_intrinsic_map,
    size_t number_of_gcps,
    KeysetContainer& keysets_gcp,
    Point3DContainer& gcps,
    hs::sfm::TrackContainer& tracks_gcp) const	//调用MultipleFlightGenerator的方法生成随机的像控点数据
  {
    scene_generator_.GeneratePoints(number_of_gcps, gcps);

    CameraViewContainer camera_view_gcp;
    if (keyset_generator_(intrinsic_params_set,
                          extrinsic_params_set_absolute,
                          images,
                          gcps,
                          image_intrinsic_map,
                          keysets_gcp,
                          tracks_gcp,
                          camera_view_gcp) != 0)	//将像控点投影至每张图片
    {
      return -1;
    }

    return 0;
  }

  //Err GenerateRelativeScene(
  //  const hs::sfm::TrackContainer& tracks,
  //  const hs::sfm::CameraViewContainer& camera_views,
  //  const ExtrinsicParamsContainer& extrinsic_params_set_absolute,
  //  const Point3DContainer& points_absolute,
  //  RMatrix& rotation_similar,
  //  Translate& translate_similar,
  //  Scalar& scale_similar,
  //  size_t& camera_id_identity,
  //  size_t& camera_id_relative,
  //  ExtrinsicParamsContainer& extrinsic_params_set_relative,
  //  Point3DContainer& points_relative) const
  //{
  //  extrinsic_params_set_relative = extrinsic_params_set_absolute;
  //  points_relative = points_absolute;

  //  return (relative_generator_(tracks,
  //                              camera_views,
  //                              extrinsic_params_set_relative,
  //                              points_relative,
  //                              rotation_similar,
  //                              translate_similar,
  //                              scale_similar,
  //                              camera_id_identity,
  //                              camera_id_relative));
  //}

  Err GenerateNoisedKeysets(const KeysetContainer& keysets_true,
                            const ImageContainer& images,
                            KeysetContainer& keysets_noised,
                            bool generate_outliers = true) const	//为图像上的同名点坐标加入随机噪声
  {
    typedef EIGEN_VECTOR(Scalar, 2) Key;
    keysets_noised = keysets_true;
    EIGEN_MATRIX(Scalar, 2, 2) covariance_key;
    covariance_key.setIdentity();
    covariance_key *= key_stddev_ * key_stddev_;
    size_t number_of_keysets = keysets_true.size();
    for (size_t i = 0; i < number_of_keysets; i++)
    {
      size_t number_of_keys = keysets_true[i].size();
      for (size_t j = 0; j < number_of_keys; j++)
      {
        Scalar random;
        hs::math::random::UniformRandomVar<Scalar, 1>::Generate(
          0, 1, random);
        if (random < outlier_ratio_ && generate_outliers)
        {
          Key min_key;
          min_key.setZero();
          Key max_key;
          max_key << Scalar(images[i].m_width),
                     Scalar(images[i].m_height);
          hs::math::random::UniformRandomVar<Scalar, 2>::Generate(
            min_key, max_key, keysets_noised[i][j]);
        }
        else
        {
          Key mean = keysets_true[i][j];
          hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
            mean, covariance_key, keysets_noised[i][j]);
        }
      }
    }

    return 0;
  }

  //Err GenerateInitialScene(
  //  const KeysetContainer& keysets_noised,
  //  const IntrinsicParamsContainer& intrinsic_params_set,
  //  const hs::sfm::TrackContainer& tracks,
  //  const hs::sfm::CameraViewContainer& camera_views,
  //  const ExtrinsicParamsContainer& extrinsic_params_set_relative,
  //  const Point3DContainer& points_relative,
  //  size_t camera_id_identity,
  //  size_t camera_id_relative,
  //  ExtrinsicParamsContainer& extrinsic_params_set_relative_estimate,
  //  Point3DContainer& points_relative_estimate,
  //  hs::sfm::ObjectIndexMap& image_extrinsic_map,
  //  hs::sfm::ObjectIndexMap& track_point_map) const
  //{
  //  extrinsic_params_set_relative_estimate.clear();
  //  extrinsic_params_set_relative_estimate.push_back(
  //    extrinsic_params_set_relative[camera_id_identity]);
  //  extrinsic_params_set_relative_estimate.push_back(
  //    extrinsic_params_set_relative[camera_id_relative]);

  //  image_extrinsic_map.Resize(extrinsic_params_set_relative.size());
  //  image_extrinsic_map[camera_id_identity] = 0;
  //  image_extrinsic_map[camera_id_relative] = 1;

  //  points_relative_estimate.clear();
  //  track_point_map.Resize(points_relative.size());
  //  const hs::sfm::CameraView& view_identity =
  //    camera_views[camera_id_identity];
  //  const hs::sfm::CameraView& view_relative =
  //    camera_views[camera_id_relative];
  //  for (size_t i = 0; i < view_identity.size(); i++)
  //  {
  //    size_t track_id = view_identity[i].first;
  //    for (size_t j = 0; j < view_relative.size(); j++)
  //    {
  //      if (track_id == view_relative[j].first)
  //      {
  //        track_point_map[track_id] = points_relative_estimate.size();
  //        points_relative_estimate.push_back(points_relative[track_id]);
  //        break;
  //      }
  //    }
  //  }

  //  return 0;
  //}

  Scalar key_stddev() const
  {
    return key_stddev_;
  }

private:
  SceneGenerator scene_generator_;
  KeysetGenerator keyset_generator_;
  //RelativeGenerator relative_generator_;
  Scalar outlier_ratio_;
  Scalar key_stddev_;
  IntrinsicParamsContainer intrinsic_params_set_true_;	//相机内参数矩阵容器
  size_t number_of_points_;
  size_t number_of_gcps_;
  size_t number_of_check_points_;
};

}
}
}

#endif
