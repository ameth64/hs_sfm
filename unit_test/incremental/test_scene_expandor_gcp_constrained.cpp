#include <iostream>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"
#include "hs_math/geometry/euler_angles.hpp"

#include "hs_sfm/synthetic/keyset_generator.hpp"
#include "hs_sfm/synthetic/relative_generator.hpp"
#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"

#include "hs_sfm/incremental/scene_expandor_gcp_constrained.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestSceneExpandorGCPConstrained
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
  typedef hs::sfm::synthetic::SceneGenerator<Scalar, ImageDimension>
          SceneGenerator;
  typedef typename SceneGenerator::IntrinsicParams IntrinsicParams;
  typedef typename IntrinsicParams::KMatrix KMatrix;
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

  typedef hs::sfm::synthetic::RelativeGenerator<Scalar> RelativeGenerator;
  typedef typename RelativeGenerator::RMatrix RMatrix;
  typedef typename RelativeGenerator::Translate Translate;

  typedef hs::sfm::incremental::SceneExpandorGCPConstrained<Scalar> Expandor;
  typedef typename Expandor::TrackPointMap TrackPointMap;
  typedef typename Expandor::ImageExtrinsicMap ImageExtrinsicMap;
  
  typedef hs::sfm::CameraFunctions<Scalar> CameraFunctions;
  typedef typename CameraFunctions::ProjectionMatrix PMatrix;

  typedef EIGEN_VECTOR(Scalar, 2) Key;

public:
  TestSceneExpandorGCPConstrained(
    Scalar focal_length_in_metre,
    size_t number_of_strips,
    size_t number_of_cameras_in_strip,
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
    Scalar camera_rotation_stddev,
    Scalar north_west_angle,
    Scalar outlier_ratio,
    Scalar key_stddev,
    const std::string& accuracy_report_path,
    const std::string& scene_data_path)
    : scene_generator_(focal_length_in_metre,
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
                       camera_plannar_stddev,
                       camera_rotation_stddev,
                       north_west_angle),
      keyset_generator_(image_width, image_height),
      outlier_ratio_(outlier_ratio),
      key_stddev_(key_stddev),
      accuracy_report_path_(accuracy_report_path),
      scene_data_path_(scene_data_path) {}

  Err Test () const
  {
    IntrinsicParamsContainer intrinsic_params_set;
    ExtrinsicParamsContainer extrinsic_params_set_absolute;
    ExtrinsicParamsContainer extrinsic_params_set_relative;
    ImageContainer images;
    Point3DContainer points_absolute;
    Point3DContainer points_relative;
    KeysetContainer keysets_true;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;
    Point3DContainer gcps;
    KeysetContainer keysets_gcp_true;
    hs::sfm::TrackContainer tracks_gcp;
    hs::sfm::CameraViewContainer camera_views_gcp;
    RMatrix rotation_similar;
    Translate translate_similar;
    Scalar scale_similar;
    size_t camera_id_identity;
    size_t camera_id_relative;
    if (TestGenerator(intrinsic_params_set,
                      extrinsic_params_set_absolute,
                      extrinsic_params_set_relative,
                      images,
                      points_absolute,
                      points_relative,
                      keysets_true,
                      tracks,
                      camera_views,
                      gcps,
                      keysets_gcp_true,
                      tracks_gcp,
                      camera_views_gcp,
                      rotation_similar,
                      translate_similar,
                      scale_similar,
                      camera_id_identity,
                      camera_id_relative) != 0)
    {
      return -1;
    }

    KeysetContainer keysets_noised;
    if (TestNoiser(keysets_true, keysets_noised, false) != 0) return -1;
    KeysetContainer keysets_gcp_noised;
    if (TestNoiser(keysets_gcp_true, keysets_gcp_noised, false) != 0) return -1;

    ExtrinsicParamsContainer extrinsic_params_set_estimate;
    Point3DContainer points_estimate;
    ImageExtrinsicMap image_extrinsic_map;
    TrackPointMap track_point_map;
    if (TestExpandor(keysets_noised,
                     intrinsic_params_set,
                     tracks,
                     camera_views,
                     extrinsic_params_set_relative,
                     points_relative,
                     gcps,
                     keysets_gcp_noised,
                     tracks_gcp,
                     camera_id_identity,
                     camera_id_relative,
                     extrinsic_params_set_estimate,
                     points_estimate,
                     image_extrinsic_map,
                     track_point_map) != 0)
    {
      return -1;
    }

    if (TestAccuracy(intrinsic_params_set,
                     extrinsic_params_set_absolute,
                     images,
                     points_absolute,
                     image_extrinsic_map,
                     track_point_map,
                     extrinsic_params_set_estimate,
                     points_estimate) != 0)
    {
      return -1;
    }

    return 0;
  }

private:
  Err TestGenerator(IntrinsicParamsContainer& intrinsic_params_set,
                    ExtrinsicParamsContainer& extrinsic_params_set_absolute,
                    ExtrinsicParamsContainer& extrinsic_params_set_relative,
                    ImageContainer& images,
                    Point3DContainer& points_absolute,
                    Point3DContainer& points_relative,
                    KeysetContainer& keysets,
                    hs::sfm::TrackContainer& tracks,
                    hs::sfm::CameraViewContainer& camera_views,
                    Point3DContainer& gcps,
                    KeysetContainer& keysets_gcp,
                    hs::sfm::TrackContainer& tracks_gcp,
                    hs::sfm::CameraViewContainer& camera_views_gcp,
                    RMatrix& rotation_similar,
                    Translate& translate_similar,
                    Scalar& scale_similar,
                    size_t& camera_id_identity,
                    size_t& camera_id_relative) const
  {
    if (scene_generator_(intrinsic_params_set,
                         extrinsic_params_set_absolute,
                         images,
                         points_absolute) != 0)
    {
      return -1;
    }

    if (keyset_generator_(intrinsic_params_set,
                          extrinsic_params_set_absolute,
                          points_absolute,
                          keysets,
                          tracks,
                          camera_views) != 0)
    {
      return -1;
    }

    size_t number_of_gcps = 7;
    scene_generator_.GenerateScenePoints(number_of_gcps, gcps);

    if (keyset_generator_(intrinsic_params_set,
                          extrinsic_params_set_absolute,
                          gcps,
                          keysets_gcp,
                          tracks_gcp,
                          camera_views_gcp) != 0)
    {
      return -1;
    }

    size_t number_of_gcps_avaiable = 0;
    for (size_t i = 0; i < number_of_gcps; i++)
    {
      if (tracks_gcp[i].size() < 4)
      {
        tracks_gcp[i].clear();
      }
      else
      {
        number_of_gcps_avaiable++;
      }
    }

    std::cout<<"number of avaiable gcps:"<<number_of_gcps_avaiable<<"\n";

    extrinsic_params_set_relative = extrinsic_params_set_absolute;
    points_relative = points_absolute;
    if (relative_generator_(tracks,
                            camera_views,
                            extrinsic_params_set_relative,
                            points_relative,
                            rotation_similar,
                            translate_similar,
                            scale_similar,
                            camera_id_identity,
                            camera_id_relative) != 0)
    {
      return -1;
    }

    return 0;
  }

  Err TestNoiser(const KeysetContainer& keysets_true,
                 KeysetContainer& keysets_noised,
                 bool has_outlier) const
  {
    keysets_noised = keysets_true;
    EIGEN_MATRIX(Scalar, 2, 2) covariance_key;
    covariance_key.setIdentity();
    covariance_key *= key_stddev_ * key_stddev_;
    Key min_key;
    min_key << -Scalar(scene_generator_.image_width()) / 2,
               -Scalar(scene_generator_.image_height()) / 2;
    Key max_key = -min_key;
    size_t number_of_keysets = keysets_true.size();
    for (size_t i = 0; i < number_of_keysets; i++)
    {
      size_t number_of_keys = keysets_true[i].size();
      for (size_t j = 0; j < number_of_keys; j++)
      {
        Scalar random;
        hs::math::random::UniformRandomVar<Scalar, 1>::Generate(
          0, 1, random);
        if (random < outlier_ratio_ && has_outlier)
        {
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

  Err TestExpandor(
    const KeysetContainer& keysets_noised,
    const IntrinsicParamsContainer& intrinsic_params_set,
    const hs::sfm::TrackContainer& tracks,
    const hs::sfm::CameraViewContainer& camera_views,
    const ExtrinsicParamsContainer& extrinsic_params_set_relative,
    const Point3DContainer& points_relative,
    const Point3DContainer& gcps,
    const KeysetContainer& keysets_noised_gcp,
    const hs::sfm::TrackContainer& tracks_gcp,
    size_t camera_id_identity,
    size_t camera_id_relative,
    ExtrinsicParamsContainer& extrinsic_params_set_estimate,
    Point3DContainer& points_estimate,
    ImageExtrinsicMap& image_extrinsic_map,
    TrackPointMap& track_point_map) const
  {
    extrinsic_params_set_estimate.clear();
    extrinsic_params_set_estimate.push_back(
      extrinsic_params_set_relative[camera_id_identity]);
    extrinsic_params_set_estimate.push_back(
      extrinsic_params_set_relative[camera_id_relative]);

    image_extrinsic_map.Resize(extrinsic_params_set_relative.size());
    image_extrinsic_map[camera_id_identity] = 0;
    image_extrinsic_map[camera_id_relative] = 1;

    points_estimate.clear();
    track_point_map.Resize(points_relative.size());
    const hs::sfm::CameraView& view_identity = camera_views[camera_id_identity];
    const hs::sfm::CameraView& view_relative = camera_views[camera_id_relative];
    for (size_t i = 0; i < view_relative.size(); i++)
    {
      size_t track_id = view_identity[i].first;
      for (size_t j = 0; j < view_relative.size(); j++)
      {
        if (track_id == view_relative[j].first)
        {
          track_point_map[track_id] = points_estimate.size();
          points_estimate.push_back(points_relative[track_id]);
          break;
        }
      }
    }

    Expandor expandor(8, key_stddev_ * 4);
    Scalar reprojection_error;
    if (expandor(keysets_noised,
                 intrinsic_params_set,
                 tracks,
                 gcps,
                 keysets_noised_gcp,
                 tracks_gcp,
                 extrinsic_params_set_estimate,
                 image_extrinsic_map,
                 points_estimate,
                 track_point_map,
                 reprojection_error) != 0)
    {
      return -1;
    }

    std::cout<<"reprojection_error:"<<reprojection_error<<"\n";
    if (reprojection_error < key_stddev_ + 1)
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }

  Err TestAccuracy(
    const IntrinsicParamsContainer& intrinsic_params_set,
    const ExtrinsicParamsContainer& extrinsic_params_set_true,
    const ImageContainer& images,
    const Point3DContainer& points_true,
    const ImageExtrinsicMap& image_extrinsic_map,
    const TrackPointMap& track_point_map,
    const ExtrinsicParamsContainer& extrinsic_params_set_estimate,
    const Point3DContainer& points_estimate) const
  {
    std::ofstream report_file(accuracy_report_path_.c_str(), std::ios::out);
    if (!report_file.is_open())
    {
      return -1;
    }

    typedef hs::math::geometry::EulerAngles<Scalar> EulerAngles;
    size_t number_of_images = images.size();
    Scalar mean_position_planar_error = Scalar(0);
    Scalar mean_position_height_error = Scalar(0);
    Scalar mean_rotation_angle0_error = Scalar(0);
    Scalar mean_rotation_angle1_error = Scalar(0);
    Scalar mean_rotation_angle2_error = Scalar(0);
    size_t number_of_extrinsic_estimate = 0;
    report_file<<"extrinsic params accuracy:\n";
    for (size_t i = 0; i < number_of_images; i++)
    {
      if (image_extrinsic_map.IsValid(i))
      {
        size_t extrinsic_id = image_extrinsic_map[i];
        const ExtrinsicParams& extrinsic_params_estimate =
          extrinsic_params_set_estimate[extrinsic_id];
        const ExtrinsicParams& extrinsic_params_true =
          extrinsic_params_set_true[i];

        RMatrix rotation_estimate =
          extrinsic_params_estimate.rotation();
        EulerAngles angles_estimate;
        angles_estimate.template FromOrthoRotMat<2, 1, -3, 1>(
          rotation_estimate);

        RMatrix rotation_true =
          extrinsic_params_true.rotation();
        EulerAngles angles_true;
        angles_true.template FromOrthoRotMat<2, 1, -3, 1>(
          rotation_true);

        Scalar rotation_angle0_error =
          std::abs(angles_true[0] - angles_estimate[0]);
        Scalar rotation_angle1_error =
          std::abs(angles_true[1] - angles_estimate[1]);
        Scalar rotation_angle2_error =
          std::abs(angles_true[2] - angles_estimate[2]);

        const Point3D& position_estimate =
          extrinsic_params_estimate.position();
        const Point3D& position_true =
          extrinsic_params_true.position();
        
        Point3D position_diff = position_estimate - position_true;
        Scalar position_planar_error = position_diff.segment(0, 2).norm();
        Scalar position_height_error = std::abs(position_diff[2]);

        mean_position_planar_error += position_planar_error;
        mean_position_height_error += position_height_error;

        mean_rotation_angle0_error += rotation_angle0_error;
        mean_rotation_angle1_error += rotation_angle1_error;
        mean_rotation_angle2_error += rotation_angle2_error;

        report_file<<i<<" "
                   <<position_estimate[0]<<" "
                   <<position_estimate[1]<<" "
                   <<position_estimate[2]<<" "
                   <<position_planar_error<<" "
                   <<position_height_error<<" "
                   <<rotation_angle0_error<<" "
                   <<rotation_angle1_error<<" "
                   <<rotation_angle2_error<<"\n";

        number_of_extrinsic_estimate++;
      }
    }
    mean_position_planar_error /= Scalar(number_of_extrinsic_estimate);
    mean_position_height_error /= Scalar(number_of_extrinsic_estimate);
    mean_rotation_angle0_error /= Scalar(number_of_extrinsic_estimate);
    mean_rotation_angle1_error /= Scalar(number_of_extrinsic_estimate);
    mean_rotation_angle2_error /= Scalar(number_of_extrinsic_estimate);
    report_file<<"mean position planar error:"
               <<mean_position_planar_error<<"\n";
    report_file<<"mean position height error:"
               <<mean_position_height_error<<"\n";
    report_file<<"mean rotation angle0 error:"
               <<mean_rotation_angle0_error<<"\n";
    report_file<<"mean rotation angle1 error:"
               <<mean_rotation_angle1_error<<"\n";
    report_file<<"mean rotation angle2 error:"
               <<mean_rotation_angle2_error<<"\n";

    size_t number_of_tracks = track_point_map.Size();
    Scalar mean_point_planar_error = Scalar(0);
    Scalar mean_point_height_error = Scalar(0);
    size_t number_of_points_estimate = 0;
    report_file<<"points accuracy:\n";
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      if (track_point_map.IsValid(i))
      {
        size_t point_id = track_point_map[i];
        const Point3D& point_estimate =
          points_estimate[point_id];
        Point3D diff = point_estimate - points_true[i];
        Scalar planar_error = diff.segment(0, 2).norm();
        Scalar height_error = std::abs(diff[2]);
        number_of_points_estimate++;
        mean_point_planar_error += planar_error;
        mean_point_height_error += height_error;

        report_file<<point_id<<" "
                   <<point_estimate[0]<<" "
                   <<point_estimate[1]<<" "
                   <<point_estimate[2]<<" "
                   <<planar_error<<" "
                   <<height_error<<"\n";
      }
    }
    mean_point_planar_error /= Scalar(number_of_points_estimate);
    mean_point_height_error /= Scalar(number_of_points_estimate);

    report_file<<"mean point planar error:"
               <<mean_point_planar_error<<"\n";
    report_file<<"mean point height error:"
               <<mean_point_height_error<<"\n";

    return 0;
  }

private:
  SceneGenerator scene_generator_;
  KeysetGenerator keyset_generator_;
  RelativeGenerator relative_generator_;
  Scalar outlier_ratio_;
  Scalar key_stddev_;
  std::string accuracy_report_path_;
  std::string scene_data_path_;
};

TEST(TestSceneExpandorGCPConstrained, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestSceneExpandorGCPConstrained<Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 3;
  size_t number_of_cameras_in_strips = 3;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 200;
  Scalar lateral_overlap_ratio = 0.7;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 100;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 5;
  Scalar north_west_angle = 60;
  Scalar outlier_ratio = 0.01;
  Scalar key_stddev = 2.0;

  Test test(focal_length_in_metre,
            number_of_strips,
            number_of_cameras_in_strips,
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
            north_west_angle,
            outlier_ratio,
            key_stddev,
            "small_data_accuracy_gcp.txt",
            "small_data_scene_gcp.xug");

  ASSERT_EQ(0, test.Test());
}

TEST(TestSceneExpandorGCPConstrained, BigDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestSceneExpandorGCPConstrained<Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 7;
  size_t number_of_cameras_in_strips = 15;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 20000;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 100;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 5;
  Scalar north_west_angle = 60;
  Scalar outlier_ratio = 0.0;
  Scalar key_stddev = 1;

  Test test(focal_length_in_metre,
            number_of_strips,
            number_of_cameras_in_strips,
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
            north_west_angle,
            outlier_ratio,
            key_stddev,
            "big_data_accuracy_gcp.txt",
            "big_data_scene_gcp.xug");

  ASSERT_EQ(0, test.Test());
}

}