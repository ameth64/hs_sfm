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

#include "hs_sfm/incremental/scene_expandor.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestSceneExpandor
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

  typedef hs::sfm::incremental::SceneExpandor<Scalar> Expandor;
  typedef typename Expandor::TrackPointMap TrackPointMap;
  typedef typename Expandor::ImageExtrinsicMap ImageExtrinsicMap;

  typedef hs::sfm::CameraFunctions<Scalar> CameraFunctions;
  typedef typename CameraFunctions::ProjectionMatrix PMatrix;

  typedef EIGEN_VECTOR(Scalar, 2) Key;

public:
  TestSceneExpandor(
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

  Err Test ()
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
                      rotation_similar,
                      translate_similar,
                      scale_similar,
                      camera_id_identity,
                      camera_id_relative) != 0)
    {
      return -1;
    }

    KeysetContainer keysets_noised;
    if (TestNoiser(keysets_true, keysets_noised) != 0) return -1;

    ExtrinsicParamsContainer extrinsic_params_set_relative_estimate;
    Point3DContainer points_relative_estimate;
    ImageExtrinsicMap image_extrinsic_map;
    TrackPointMap track_point_map;
    if (TestExpandor(keysets_noised,
                    intrinsic_params_set,
                    tracks,
                    camera_views,
                    extrinsic_params_set_relative,
                    points_relative,
                    camera_id_identity,
                    camera_id_relative,
                    extrinsic_params_set_relative_estimate,
                    points_relative_estimate,
                    image_extrinsic_map,
                    track_point_map) != 0)
    {
      return -1;
    }

    if (TestAccuracy(intrinsic_params_set,
                     extrinsic_params_set_absolute,
                     images,
                     points_absolute,
                     rotation_similar,
                     translate_similar,
                     scale_similar,
                     image_extrinsic_map,
                     track_point_map,
                     extrinsic_params_set_relative_estimate,
                     points_relative_estimate) != 0)
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

    std::ofstream scene_file(scene_data_path_.c_str(), std::ios::out);
    size_t number_of_cameras = extrinsic_params_set_absolute.size();
    scene_file<<"camera data:\n";
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      const ExtrinsicParams& extrinsic_params =
        extrinsic_params_set_absolute[i];
      scene_file<<i<<" "
                <<extrinsic_params.position()[0]<<" "
                <<extrinsic_params.position()[1]<<" "
                <<extrinsic_params.position()[2]<<"\n";
    }
    size_t number_of_points = points_absolute.size();
    scene_file<<"point data:\n";
    for (size_t i = 0; i < number_of_points; i++)
    {
      scene_file<<i<<" "
                <<points_absolute[i][0]<<" "
                <<points_absolute[i][1]<<" "
                <<points_absolute[i][2]<<"\n";
    }

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
                 KeysetContainer& keysets_noised) const
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
        if (random < outlier_ratio_)
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
    size_t camera_id_identity,
    size_t camera_id_relative,
    ExtrinsicParamsContainer& extrinsic_params_set_relative_estimate,
    Point3DContainer& points_relative_estimate,
    ImageExtrinsicMap& image_extrinsic_map,
    TrackPointMap& track_point_map) const
  {
    extrinsic_params_set_relative_estimate.clear();
    extrinsic_params_set_relative_estimate.push_back(
      extrinsic_params_set_relative[camera_id_identity]);
    extrinsic_params_set_relative_estimate.push_back(
      extrinsic_params_set_relative[camera_id_relative]);

    image_extrinsic_map.Resize(extrinsic_params_set_relative.size());
    image_extrinsic_map[camera_id_identity] = 0;
    image_extrinsic_map[camera_id_relative] = 1;

    points_relative_estimate.clear();
    track_point_map.Resize(points_relative.size());
    const hs::sfm::CameraView& view_identity = camera_views[camera_id_identity];
    const hs::sfm::CameraView& view_relative = camera_views[camera_id_relative];
    for (size_t i = 0; i < view_identity.size(); i++)
    {
      size_t track_id = view_identity[i].first;
      for (size_t j = 0; j < view_relative.size(); j++)
      {
        if (track_id == view_relative[j].first)
        {
          track_point_map[track_id] = points_relative_estimate.size();
          points_relative_estimate.push_back(points_relative[track_id]);
          break;
        }
      }
    }

    Expandor expandor(8, key_stddev_ * 4);
    Scalar reprojection_error;
    if (expandor(keysets_noised, intrinsic_params_set, tracks,
                 extrinsic_params_set_relative_estimate,
                 image_extrinsic_map,
                 points_relative_estimate,
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
    const ExtrinsicParamsContainer& extrinsic_params_set_absolute,
    const ImageContainer& images,
    const Point3DContainer& points_absolute,
    const RMatrix& rotation_similar,
    const Translate translate_similar,
    Scalar scale_similar,
    const ImageExtrinsicMap& image_extrinsic_map,
    const TrackPointMap& track_point_map,
    const ExtrinsicParamsContainer& extrinsic_params_set_relative_estimate,
    const Point3DContainer& points_relative_estimate) const
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
    for (size_t i = 0 ; i < number_of_images; i++)
    {
      if (image_extrinsic_map.IsValid(i))
      {
        size_t extrinsic_id = image_extrinsic_map[i];
        const ExtrinsicParams& extrinsic_params_relative_estimate =
          extrinsic_params_set_relative_estimate[extrinsic_id];
        const ExtrinsicParams& extrinsic_params_absolute =
          extrinsic_params_set_absolute[i];

        RMatrix rotation_relative_estimate =
          extrinsic_params_relative_estimate.rotation();
        RMatrix rotation_absolute_estimate =
          rotation_relative_estimate * rotation_similar.transpose();
        EulerAngles angles_absolute_estimate;
        angles_absolute_estimate.template FromOrthoRotMat<2, 1, -3, 1>(
          rotation_absolute_estimate);

        RMatrix rotation_absolute = extrinsic_params_absolute.rotation();
        EulerAngles angles_absolute;
        angles_absolute.template FromOrthoRotMat<2, 1, -3, 1>(
          rotation_absolute);

        Scalar rotation_angle0_error =
          std::abs(angles_absolute[0] - angles_absolute_estimate[0]);
        Scalar rotation_angle1_error =
          std::abs(angles_absolute[1] - angles_absolute_estimate[1]);
        Scalar rotation_angle2_error =
          std::abs(angles_absolute[2] - angles_absolute_estimate[2]);

        const Point3D& position_relative_estimate =
          extrinsic_params_relative_estimate.position();
        Point3D position_absolute_estimate =
          scale_similar * rotation_similar * position_relative_estimate +
          translate_similar;
        const Point3D& position_absolute =
          extrinsic_params_absolute.position();

        Point3D position_diff = position_absolute_estimate - position_absolute;
        Scalar position_planar_error = position_diff.segment(0, 2).norm();
        Scalar position_height_error = std::abs(position_diff[2]);

        mean_position_planar_error += position_planar_error;
        mean_position_height_error += position_height_error;

        mean_rotation_angle0_error += rotation_angle0_error;
        mean_rotation_angle1_error += rotation_angle1_error;
        mean_rotation_angle2_error += rotation_angle2_error;

        report_file<<i<<" "
                   <<position_absolute_estimate[0]<<" "
                   <<position_absolute_estimate[1]<<" "
                   <<position_absolute_estimate[2]<<" "
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
        const Point3D& point_relative_estimate =
          points_relative_estimate[point_id];
        Point3D point_absolute_estimate =
          scale_similar * rotation_similar * point_relative_estimate +
          translate_similar;
        Point3D diff = point_absolute_estimate - points_absolute[i];
        Scalar planar_error = diff.segment(0, 2).norm();
        Scalar height_error = std::abs(diff[2]);
        number_of_points_estimate++;
        mean_point_planar_error += planar_error;
        mean_point_height_error += height_error;

        report_file<<point_id<<" "
                   <<point_absolute_estimate[0]<<" "
                   <<point_absolute_estimate[1]<<" "
                   <<point_absolute_estimate[2]<<" "
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

TEST(TestSceneExpandor, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestSceneExpandor<Scalar, ImageDimension> Test;

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
            "small_data_accuracy.txt",
            "small_data_scene.xug");

  ASSERT_EQ(0, test.Test());
}

TEST(TestSceneExpandor, BigDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestSceneExpandor<Scalar, ImageDimension> Test;

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
            "big_data_accuracy.txt",
            "big_data_scene.xug");

  ASSERT_EQ(0, test.Test());
}

}
