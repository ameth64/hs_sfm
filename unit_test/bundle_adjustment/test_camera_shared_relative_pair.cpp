#include <iostream>
#include <fstream>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/geometry/euler_angles.hpp"
#include "hs_sfm/sfm_utility/similar_transform_estimator.hpp"
#include "hs_sfm/sfm_utility/projective_functions.hpp"
#include "hs_sfm/synthetic/relative_pair_generator.hpp"
#include "hs_sfm/synthetic/multiple_camera_keyset_generator.hpp"
#include "hs_sfm/essential/ematrix_5_points_calculator.hpp"
#include "hs_sfm/essential/ematrix_extrinsic_params_points_calculator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_ceres_optimizor.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestCameraSharedRelativePair
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

  typedef hs::sfm::synthetic::RelativePairGenerator<Scalar, ImageDimension>
          RelativePairGenerator;
  typedef typename RelativePairGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename RelativePairGenerator::Image Image;
  typedef typename RelativePairGenerator::ImageContainer ImageContainer;
  typedef typename RelativePairGenerator::Point3D Point3D;
  typedef typename RelativePairGenerator::Point3DContainer Point3DContainer;

  typedef hs::sfm::synthetic::MultipleCameraKeysetGenerator<Scalar,
                                                            ImageDimension>
          KeysetGenerator;
  typedef typename KeysetGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename KeysetGenerator::IntrinsicParams IntrinsicParams;
  typedef typename IntrinsicParams::KMatrix KMatrix;
  typedef typename KeysetGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

  typedef hs::sfm::essential::EMatrix5PointsCalculator<Scalar>
          EMatrixCalculator;
  typedef typename EMatrixCalculator::EMatrix EMatrix;
  typedef typename EMatrixCalculator::EMatrixHypotheses EMatrixHypotheses;
  typedef typename EMatrixCalculator::HKey HKey;
  typedef typename EMatrixCalculator::HKeyPair HKeyPair;
  typedef typename EMatrixCalculator::HKeyPairContainer HKeyPairContainer;
  typedef hs::sfm::essential::EMatrixExtrinsicParamsPointsCalculator<Scalar>
          ExtrinsicParamPointsCalculator;

  typedef hs::sfm::SimilarTransformEstimator<Scalar> SimilarEstimator;
  typedef typename SimilarEstimator::Rotation Rotation;
  typedef typename SimilarEstimator::Translate Translate;

  typedef hs::math::geometry::EulerAngles<Scalar> EulerAngles;

  typedef hs::sfm::ba::CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::FeatureMap FeatureMap;
  typedef typename VectorFunction::FeatureMapContainer FeatureMapContainer;
  typedef typename VectorFunction::ImageCameraMap ImageCameraMap;
  typedef typename VectorFunction::Image BAImage;
  typedef typename VectorFunction::ImageContainer BAImageContainer;
  typedef typename VectorFunction::Camera BACamera;
  typedef typename VectorFunction::CameraContainer BACameraContainer;
  typedef hs::sfm::ba::CameraSharedLevenbergMarquardtOptimizor<VectorFunction>
          Optimizor;
  typedef typename Optimizor::YCovarianceInverse YCovarianceInverse;
  typedef hs::sfm::ba::CameraSharedCeresOptimizor<VectorFunction>
          CeresOptimizor;

private:
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_MATRIX(Scalar, 2, 2) Matrix22;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

public:
  TestCameraSharedRelativePair(
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    const IntrinsicParams& intrinsic_params_identity,
    const IntrinsicParams& intrinsic_params_relative,
    Scalar lateral_overlap_ratio,
    Scalar longitudinal_overlap_ratio,
    Scalar scene_max_height,
    Scalar camera_height_stddev,
    Scalar camera_planar_stddev,
    Scalar camera_rotation_stddev,
    size_t number_of_points,
    bool is_uniform_camera,
    Scalar key_stddev)
    : relative_pair_generator_(ground_resolution,
                               image_width,
                               image_height,
                               pixel_size,
                               intrinsic_params_identity,
                               is_uniform_camera ?
                               intrinsic_params_identity :
                               intrinsic_params_relative,
                               lateral_overlap_ratio,
                               longitudinal_overlap_ratio,
                               scene_max_height,
                               camera_height_stddev,
                               camera_planar_stddev,
                               camera_rotation_stddev),
      keyset_generator_(),
      number_of_points_(number_of_points),
      is_uniform_camera_(is_uniform_camera),
      key_stddev_(key_stddev) {}

  Err Test() const
  {
    ExtrinsicParams extrinsic_params_identity_true,
                    extrinsic_params_relative_true;
    Point3DContainer points_abs;
    KeysetContainer keysets;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;
    if (GenerateTestPair(extrinsic_params_identity_true,
                         extrinsic_params_relative_true,
                         points_abs,
                         keysets,
                         tracks,
                         camera_views) != 0)
    {
      return -1;
    }

    const IntrinsicParams& intrinsic_params_identity_true =
      relative_pair_generator_.intrinsic_params_identity();
    const IntrinsicParams& intrinsic_params_relative_true =
      is_uniform_camera_ ?
      relative_pair_generator_.intrinsic_params_identity() :
      relative_pair_generator_.intrinsic_params_relative();

    IntrinsicParams intrinsic_params_identity_essential;
    intrinsic_params_identity_essential.set_focal_length(
      intrinsic_params_identity_true.focal_length() + 10);
    IntrinsicParams intrinsic_params_relative_essential;
    if (is_uniform_camera_)
    {
      intrinsic_params_relative_essential.set_focal_length(
        intrinsic_params_identity_true.focal_length() + 10);
    }
    else
    {
      intrinsic_params_relative_essential.set_focal_length(
        intrinsic_params_relative_true.focal_length() + 10);
    }

    KeysetContainer keysets_noised;
    if (GenerateNoisedKeysets(keysets, keysets_noised) != 0)
    {
      return -1;
    }

    HKeyPairContainer hkey_pairs;
    if (GenerateHomogenerousKeyPairs(keysets_noised,
                                     intrinsic_params_identity_essential,
                                     intrinsic_params_relative_essential,
                                     hkey_pairs) != 0)
    {
      return -1;
    }

    ExtrinsicParams extrinsic_params_relative_essential;
    Point3DContainer points_essential;
    if (ComputeEssentialExtrinsicParamsPoints(
          hkey_pairs,
          extrinsic_params_relative_essential,
          points_essential) != 0)
    {
      return -1;
    }

    Rotation rotation_similar;
    Translate translate_similar;
    Scalar scale_similar;
    ExtrinsicParams extrinsic_params_identity_essential;
    Matrix33 ortho_matrix_identity = Matrix33::Identity();
    Vector3 zero = Vector3::Zero();
    extrinsic_params_identity_essential.rotation() = ortho_matrix_identity;
    extrinsic_params_identity_essential.position() = zero;

    std::cout<<"Essential Error:\n";
    if (TestError(points_essential, points_abs,
                  intrinsic_params_identity_essential,
                  intrinsic_params_relative_essential,
                  extrinsic_params_identity_essential,
                  extrinsic_params_relative_essential,
                  intrinsic_params_identity_true,
                  intrinsic_params_relative_true,
                  extrinsic_params_identity_true,
                  extrinsic_params_relative_true,
                  keysets_noised,
                  rotation_similar,
                  translate_similar,
                  scale_similar) != 0)
    {
      return -1;
    }

    //旋转向量为0时（即没有旋转），bundle adjustment时无法计算jacobian矩阵。
    //因此需改变初始的相机朝向
    Rotation rotation_extra;
    rotation_extra[0] = Scalar(0);
    rotation_extra[1] = Scalar(0);
    rotation_extra[2] = Scalar(3.141592653) / Scalar(180) * Scalar(45);

    extrinsic_params_identity_essential.rotation() =
      extrinsic_params_identity_essential.rotation() * rotation_extra.Inverse();
    extrinsic_params_identity_essential.position() =
      rotation_extra * extrinsic_params_identity_essential.position();
    extrinsic_params_relative_essential.rotation() =
      extrinsic_params_relative_essential.rotation() * rotation_extra.Inverse();
    extrinsic_params_relative_essential.position() =
      rotation_extra * extrinsic_params_relative_essential.position();
    for (size_t i = 0; i < points_essential.size(); i++)
    {
      Point3D point = points_essential[i];
      points_essential[i] = rotation_extra * point;
    }

    ExtrinsicParams extrinsic_params_identity_bundle,
                    extrinsic_params_relative_bundle;
    IntrinsicParams intrinsic_params_identity_bundle,
                    intrinsic_params_relative_bundle;
    Point3DContainer points_bundle;
    hs::sfm::ba::FixMask fix_mask;
    //fix_mask.set(hs::sfm::ba::FIX_POINTS);
    //fix_mask.set(hs::sfm::ba::FIX_IMAGES);
    if (BundleAdjustment(extrinsic_params_identity_essential,
                         extrinsic_params_relative_essential,
                         intrinsic_params_identity_essential,
                         intrinsic_params_relative_essential,
                         points_essential,
                         keysets_noised,
                         fix_mask,
                         extrinsic_params_identity_bundle,
                         extrinsic_params_relative_bundle,
                         intrinsic_params_identity_bundle,
                         intrinsic_params_relative_bundle,
                         points_bundle) != 0)
    {
      return -1;
    }
    std::cout<<"Bundle Error:\n";
    if (TestError(points_bundle, points_abs,
                  intrinsic_params_identity_bundle,
                  intrinsic_params_relative_bundle,
                  extrinsic_params_identity_bundle,
                  extrinsic_params_relative_bundle,
                  intrinsic_params_identity_true,
                  intrinsic_params_relative_true,
                  extrinsic_params_identity_true,
                  extrinsic_params_relative_true,
                  keysets_noised,
                  rotation_similar,
                  translate_similar,
                  scale_similar) != 0)
    {
      return -1;
    }

    ExtrinsicParams extrinsic_params_identity_ceres,
                    extrinsic_params_relative_ceres;
    IntrinsicParams intrinsic_params_identity_ceres,
                    intrinsic_params_relative_ceres;
    Point3DContainer points_ceres;
    hs::sfm::ba::FixMask fix_mask_ceres;
    //fix_mask_ceres.set(hs::sfm::ba::FIX_POINTS);
    //fix_mask_ceres.set(hs::sfm::ba::FIX_IMAGES);
    if (BundleAdjustment(extrinsic_params_identity_essential,
                         extrinsic_params_relative_essential,
                         intrinsic_params_identity_essential,
                         intrinsic_params_relative_essential,
                         points_essential,
                         keysets_noised,
                         fix_mask_ceres,
                         extrinsic_params_identity_ceres,
                         extrinsic_params_relative_ceres,
                         intrinsic_params_identity_ceres,
                         intrinsic_params_relative_ceres,
                         points_ceres, true) != 0)
    {
      return -1;
    }
    std::cout<<"Ceres Error:\n";
    if (TestError(points_ceres, points_abs,
                  intrinsic_params_identity_ceres,
                  intrinsic_params_relative_ceres,
                  extrinsic_params_identity_ceres,
                  extrinsic_params_relative_ceres,
                  intrinsic_params_identity_true,
                  intrinsic_params_relative_true,
                  extrinsic_params_identity_true,
                  extrinsic_params_relative_true,
                  keysets_noised,
                  rotation_similar,
                  translate_similar,
                  scale_similar) != 0)
    {
      return -1;
    }
    return 0;
  }

private:
  Err GenerateTestPair(ExtrinsicParams& extrinsic_params_identity,
                       ExtrinsicParams& extrinsic_params_relative,
                       Point3DContainer& points,
                       KeysetContainer& keysets,
                       hs::sfm::TrackContainer& tracks,
                       hs::sfm::CameraViewContainer& camera_views) const
  {
    Image image_identity, image_relative;
    if (relative_pair_generator_.GenerateExtrinsicParams(
          extrinsic_params_identity,
          extrinsic_params_relative,
          image_identity,
          image_relative) != 0)
    {
      return -1;
    }

    Point3DContainer points_loose;
    if (relative_pair_generator_.GeneratePoints(number_of_points_,
                                                points_loose) != 0)
    {
      return -1;
    }

    IntrinsicParamsContainer intrinsic_params_set(2);
    intrinsic_params_set[0] =
      relative_pair_generator_.intrinsic_params_identity();
    if (is_uniform_camera_)
    {
      intrinsic_params_set[1] =
        relative_pair_generator_.intrinsic_params_identity();
    }
    else
    {
      intrinsic_params_set[1] =
        relative_pair_generator_.intrinsic_params_relative();
    }

    ExtrinsicParamsContainer extrinsic_params_set(2);
    extrinsic_params_set[0] = extrinsic_params_identity;
    extrinsic_params_set[1] = extrinsic_params_relative;

    ImageContainer images(2);
    images[0] = image_identity;
    images[1] = image_relative;

    std::vector<size_t> image_intrinsic_map(2);
    image_intrinsic_map[0] = 0;
    if (is_uniform_camera_)
    {
      image_intrinsic_map[1] = 0;
    }
    else
    {
      image_intrinsic_map[1] = 1;
    }

    KeysetContainer keysets_loose;
    hs::sfm::TrackContainer tracks_loose;
    hs::sfm::CameraViewContainer camera_views_loose;
    if (keyset_generator_(intrinsic_params_set,
                          extrinsic_params_set,
                          images,
                          points_loose,
                          image_intrinsic_map,
                          keysets_loose,
                          tracks_loose,
                          camera_views_loose) != 0)
    {
      return -1;
    }

    auto itr_track_loose = tracks_loose.begin();
    auto itr_track_loose_end = tracks_loose.end();
    size_t number_of_available_points = 0;
    for (; itr_track_loose != itr_track_loose_end; ++itr_track_loose)
    {
      if (itr_track_loose->size() == 2)
      {
        number_of_available_points++;
      }
    }
    points.clear();
    tracks.clear();
    keysets.resize(2);
    keysets[0].resize(number_of_available_points);
    keysets[1].resize(number_of_available_points);
    camera_views.resize(2);
    camera_views[0].resize(number_of_available_points);
    camera_views[1].resize(number_of_available_points);
    size_t available_point_id = 0;
    auto itr_point_loose = points_loose.begin();
    itr_track_loose = tracks_loose.begin();
    points.resize(number_of_available_points);
    tracks.resize(number_of_available_points);
    for (; itr_track_loose != itr_track_loose_end;
         ++itr_track_loose, ++itr_point_loose)
    {
      if (itr_track_loose->size() == 2)
      {
        size_t key_id_0 = (*itr_track_loose)[0].second;
        size_t key_id_1 = (*itr_track_loose)[1].second;
        points[available_point_id] = *itr_point_loose;
        keysets[0][available_point_id] = keysets_loose[0][key_id_0];
        keysets[1][available_point_id] = keysets_loose[1][key_id_1];
        tracks[available_point_id].push_back(std::make_pair(
                                               0, available_point_id));
        tracks[available_point_id].push_back(std::make_pair(
                                               1, available_point_id));
        camera_views[0][available_point_id] =
          std::make_pair(available_point_id, key_id_0);
        camera_views[1][available_point_id] =
          std::make_pair(available_point_id, key_id_1);
        available_point_id++;
      }
    }

    return 0;
  }

  Err GenerateNoisedKeysets(const KeysetContainer& keysets,
                            KeysetContainer& keysets_noised) const
  {
    Matrix22 key_covariance = Matrix22::Identity();
    key_covariance *= key_stddev_ * key_stddev_;

    keysets_noised.resize(2);

    for (size_t i = 0; i < 2; i++)
    {
      size_t number_of_keys = keysets[i].size();
      keysets_noised[i].resize(number_of_keys);
      for (size_t j = 0; j < number_of_keys; j++)
      {
        hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
          keysets[i][j], key_covariance, keysets_noised[i][j]);
      }
    }

    return 0;
  }

  Err GenerateHomogenerousKeyPairs(
    const KeysetContainer& keysets,
    const IntrinsicParams& intrinsic_params_identity,
    const IntrinsicParams& intrinsic_params_relative,
    HKeyPairContainer& hkey_pairs) const
  {
    if (keysets.size() != 2)
    {
      return -1;
    }
    size_t number_of_keys = keysets[0].size();
    if (keysets[1].size() != number_of_keys)
    {
      return -1;
    }

    KMatrix K0_inverse = intrinsic_params_identity.GetKMatrix().inverse();
    KMatrix K1_inverse = intrinsic_params_relative.GetKMatrix().inverse();

    hkey_pairs.resize(number_of_keys);
    for (size_t i = 0; i < number_of_keys; i++)
    {
      HKeyPair hkey_pair;
      hkey_pair.first.template segment<2>(0) = keysets[0][i];
      hkey_pair.second.template segment<2>(0) = keysets[1][i];
      hkey_pair.first(2) = Scalar(1);
      hkey_pair.second(2) = Scalar(1);
      hkey_pair.first = K0_inverse * hkey_pair.first;
      hkey_pair.second = K1_inverse * hkey_pair.second;
      hkey_pairs[i] = hkey_pair;
    }

    return 0;
  }

  Err ComputeEssentialExtrinsicParamsPoints(
    const HKeyPairContainer& hkey_pairs,
    ExtrinsicParams& extrinsic_params,
    Point3DContainer& points) const
  {
    EMatrix ematrix;
    EMatrixCalculator ematrix_calculator;
    if (ematrix_calculator(hkey_pairs, ematrix) != 0)
    {
      return -1;
    }

    ExtrinsicParamPointsCalculator extrinsic_params_points_calculator;
    if (extrinsic_params_points_calculator(ematrix,
                                           hkey_pairs,
                                           extrinsic_params,
                                           points) != 0)
    {
      return -1;
    }

    return 0;
  }

  Err TestError(const Point3DContainer& points_estimate,
                const Point3DContainer& points_true,
                const IntrinsicParams& intrinsic_params_identity_estimate,
                const IntrinsicParams& intrinsic_params_relative_estimate,
                const ExtrinsicParams& extrinsic_params_identity_estimate,
                const ExtrinsicParams& extrinsic_params_relative_estimate,
                const IntrinsicParams& intrinsic_params_identity_true,
                const IntrinsicParams& intrinsic_params_relative_true,
                const ExtrinsicParams& extrinsic_params_identity_true,
                const ExtrinsicParams& extrinsic_params_relative_true,
                const KeysetContainer& keysets,
                Rotation& rotation_similar,
                Translate& translate_similar,
                Scalar& scale_similar) const
  {
    Scalar mean_reprojection_error = Scalar(0);
    for (size_t i = 0; i < points_estimate.size(); i++)
    {
      const Point3D& point = points_estimate[i];
      Vector2 key_identity =
        hs::sfm::ProjectiveFunctions<Scalar>::WorldPointProjectToImageKey(
          intrinsic_params_identity_estimate,
          extrinsic_params_identity_estimate,
          point);
      Vector2 key_relative =
        hs::sfm::ProjectiveFunctions<Scalar>::WorldPointProjectToImageKey(
          intrinsic_params_relative_estimate,
          extrinsic_params_relative_estimate,
          point);

      Vector2 diff_identity = key_identity - keysets[0][i];
      Vector2 diff_relative = key_relative - keysets[1][i];
      //mean_reprojection_error += diff_identity.squaredNorm() +
      //                           diff_relative.squaredNorm();
      mean_reprojection_error += diff_identity.norm() +
                                 diff_relative.norm();
    }
    mean_reprojection_error /= Scalar(points_estimate.size() * 2);
    std::cout<<"mean_reprojection_error:"<<mean_reprojection_error<<"\n";

    SimilarEstimator similar_estimator;
    if (similar_estimator(points_estimate,
                          points_true,
                          rotation_similar,
                          translate_similar,
                          scale_similar) != 0)
    {
      return -1;
    }

    Scalar point_planar_mean = Scalar(0);
    Scalar point_height_mean = Scalar(0);
    for (size_t i = 0; i < points_true.size(); i++)
    {
      Point3D point_estimate = points_estimate[i];
      point_estimate =
        scale_similar * (rotation_similar * point_estimate) +
        translate_similar;
      Point3D diff = points_true[i] - point_estimate;
      point_planar_mean += diff.segment(0, 2).norm();
      point_height_mean += std::abs(diff[2]);
    }
    point_planar_mean /= Scalar(points_true.size());
    point_height_mean /= Scalar(points_true.size());

    std::cout<<"point_planar_mean:"<<point_planar_mean<<"\n";
    std::cout<<"point_height_mean:"<<point_height_mean<<"\n";

    Matrix33 ortho_matrix_similar = rotation_similar;

    Matrix33 rotation_estimate_identity =
      extrinsic_params_identity_estimate.rotation();
    rotation_estimate_identity =
      rotation_estimate_identity * ortho_matrix_similar.transpose();
    EulerAngles angles_estimate_identity;
    angles_estimate_identity.template FromOrthoRotMat<2, 1, -3, 1>(
      rotation_estimate_identity);
    Point3D position_estimate_identity =
      extrinsic_params_identity_estimate.position();
    position_estimate_identity =
      scale_similar * (rotation_similar * position_estimate_identity) +
      translate_similar;

    Matrix33 rotation_estimate_relative =
      extrinsic_params_relative_estimate.rotation();
    rotation_estimate_relative =
      rotation_estimate_relative * ortho_matrix_similar.transpose();
    EulerAngles angles_estimate_relative;
    angles_estimate_relative.template FromOrthoRotMat<2, 1, -3, 1>(
      rotation_estimate_relative);
    Point3D position_estimate_relative =
      extrinsic_params_relative_estimate.position();
    position_estimate_relative =
      scale_similar * (rotation_similar * position_estimate_relative) +
      translate_similar;

    Matrix33 rotation_true_identity = extrinsic_params_identity_true.rotation();
    EulerAngles angles_true_identity;
    angles_true_identity.template FromOrthoRotMat<2, 1, -3, 1>(
      rotation_true_identity);
    Point3D position_true_identity = extrinsic_params_identity_true.position();

    Matrix33 rotation_true_relative = extrinsic_params_relative_true.rotation();
    EulerAngles angles_true_relative;
    angles_true_relative.template FromOrthoRotMat<2, 1, -3, 1>(
      rotation_true_relative);
    Point3D position_true_relative = extrinsic_params_relative_true.position();

    Vector3 rotation_diff_identity;
    rotation_diff_identity<<angles_estimate_identity[0] -
                            angles_true_identity[0],
                            angles_estimate_identity[1] -
                            angles_true_identity[1],
                            angles_estimate_identity[2] -
                            angles_true_identity[2];
    Vector3 rotation_diff_relative;
    rotation_diff_relative<<angles_estimate_relative[0] -
                            angles_true_relative[0],
                            angles_estimate_relative[1] -
                            angles_true_relative[1],
                            angles_estimate_relative[2] -
                            angles_true_relative[2];
    rotation_diff_identity *= Scalar(180) / Scalar(M_PI);
    rotation_diff_relative *= Scalar(180) / Scalar(M_PI);

    Vector3 position_diff_identity = position_estimate_identity -
                                     position_true_identity;
    Vector3 position_diff_relative = position_estimate_relative -
                                     position_true_relative;

    std::cout<<"rotation_diff_identity:\n"<<rotation_diff_identity<<"\n";
    std::cout<<"rotation_diff_relative:\n"<<rotation_diff_relative<<"\n";
    std::cout<<"position_diff_identity:\n"<<position_diff_identity<<"\n";
    std::cout<<"position_diff_relative:\n"<<position_diff_relative<<"\n";

    std::cout<<"k1_identity_true:\n"
             <<intrinsic_params_identity_true.k1()<<"\n";
    std::cout<<"k1_identity_estimate:\n"
             <<intrinsic_params_identity_estimate.k1()<<"\n";
    std::cout<<"k2_identity_true:\n"
             <<intrinsic_params_identity_true.k2()<<"\n";
    std::cout<<"k2_identity_estimate:\n"
             <<intrinsic_params_identity_estimate.k2()<<"\n";
    std::cout<<"k3_identity_true:\n"
             <<intrinsic_params_identity_true.k3()<<"\n";
    std::cout<<"k3_identity_estimate:\n"
             <<intrinsic_params_identity_estimate.k3()<<"\n";
    std::cout<<"d1_identity_true:\n"
             <<intrinsic_params_identity_true.d1()<<"\n";
    std::cout<<"d1_identity_estimate:\n"
             <<intrinsic_params_identity_estimate.d1()<<"\n";
    std::cout<<"d2_identity_true:\n"
             <<intrinsic_params_identity_true.d2()<<"\n";
    std::cout<<"d2_identity_estimate:\n"
             <<intrinsic_params_identity_estimate.d2()<<"\n";
    std::cout<<"focal_length_identity_true:\n"
             <<intrinsic_params_identity_true.focal_length()<<"\n";
    std::cout<<"focal_length_identity_estimate:\n"
             <<intrinsic_params_identity_estimate.focal_length()<<"\n";
    std::cout<<"skew_identity_true:\n"
             <<intrinsic_params_identity_true.skew()<<"\n";
    std::cout<<"skew_identity_estimate:\n"
             <<intrinsic_params_identity_estimate.skew()<<"\n";
    std::cout<<"principal_point_x_identity_true:\n"
             <<intrinsic_params_identity_true.principal_point_x()<<"\n";
    std::cout<<"principal_point_x_identity_estimate:\n"
             <<intrinsic_params_identity_estimate.principal_point_x()<<"\n";
    std::cout<<"principal_point_y_identity_true:\n"
             <<intrinsic_params_identity_true.principal_point_y()<<"\n";
    std::cout<<"principal_point_y_identity_estimate:\n"
             <<intrinsic_params_identity_estimate.principal_point_y()<<"\n";
    std::cout<<"pixel_ratio_identity_true:\n"
             <<intrinsic_params_identity_true.pixel_ratio()<<"\n";
    std::cout<<"pixel_ratio_identity_estimate:\n"
             <<intrinsic_params_identity_estimate.pixel_ratio()<<"\n";

    if (!is_uniform_camera_)
    {
      std::cout<<"k1_relative_true:\n"
               <<intrinsic_params_relative_true.k1()<<"\n";
      std::cout<<"k1_relative_estimate:\n"
               <<intrinsic_params_relative_estimate.k1()<<"\n";
      std::cout<<"k2_relative_true:\n"
               <<intrinsic_params_relative_true.k2()<<"\n";
      std::cout<<"k2_relative_estimate:\n"
               <<intrinsic_params_relative_estimate.k2()<<"\n";
      std::cout<<"k3_relative_true:\n"
               <<intrinsic_params_relative_true.k3()<<"\n";
      std::cout<<"k3_relative_estimate:\n"
               <<intrinsic_params_relative_estimate.k3()<<"\n";
      std::cout<<"d1_relative_true:\n"
               <<intrinsic_params_relative_true.d1()<<"\n";
      std::cout<<"d1_relative_estimate:\n"
               <<intrinsic_params_relative_estimate.d1()<<"\n";
      std::cout<<"d2_relative_true:\n"
               <<intrinsic_params_relative_true.d2()<<"\n";
      std::cout<<"d2_relative_estimate:\n"
               <<intrinsic_params_relative_estimate.d2()<<"\n";
      std::cout<<"focal_length_relative_true:\n"
               <<intrinsic_params_relative_true.focal_length()<<"\n";
      std::cout<<"focal_length_relative_estimate:\n"
               <<intrinsic_params_relative_estimate.focal_length()<<"\n";
      std::cout<<"skew_relative_true:\n"
               <<intrinsic_params_relative_true.skew()<<"\n";
      std::cout<<"skew_relative_estimate:\n"
               <<intrinsic_params_relative_estimate.skew()<<"\n";
      std::cout<<"principal_point_x_relative_true:\n"
               <<intrinsic_params_relative_true.principal_point_x()<<"\n";
      std::cout<<"principal_point_x_relative_estimate:\n"
               <<intrinsic_params_relative_estimate.principal_point_x()<<"\n";
      std::cout<<"principal_point_y_relative_true:\n"
               <<intrinsic_params_relative_true.principal_point_y()<<"\n";
      std::cout<<"principal_point_y_relative_estimate:\n"
               <<intrinsic_params_relative_estimate.principal_point_y()<<"\n";
      std::cout<<"pixel_ratio_relative_true:\n"
               <<intrinsic_params_relative_true.pixel_ratio()<<"\n";
      std::cout<<"pixel_ratio_relative_estimate:\n"
               <<intrinsic_params_relative_estimate.pixel_ratio()<<"\n";
    }

    return 0;
  }

  Err BundleAdjustment(const ExtrinsicParams& extrinsic_params_identity_initial,
                       const ExtrinsicParams& extrinsic_params_relative_initial,
                       const IntrinsicParams& intrinsic_params_identity_initial,
                       const IntrinsicParams& intrinsic_params_relative_initial,
                       const Point3DContainer& points_initial,
                       const KeysetContainer& keysets,
                       const hs::sfm::ba::FixMask& fix_mask,
                       ExtrinsicParams& extrinsic_params_identity_bundle,
                       ExtrinsicParams& extrinsic_params_relative_bundle,
                       IntrinsicParams& intrinsic_params_identity_bundle,
                       IntrinsicParams& intrinsic_params_relative_bundle,
                       Point3DContainer& points_bundle,
                       bool use_ceres = false) const
  {
    using namespace hs::sfm::ba;
    size_t number_of_points = points_initial.size();
    size_t number_of_keys = number_of_points * 2;
    if (keysets.size() != 2 ||
        keysets[0].size() != number_of_points ||
        keysets[1].size() != number_of_points) return -1;
    VectorFunction vector_function;
    vector_function.set_number_of_images(2);
    vector_function.set_number_of_points(Index(number_of_points));
    vector_function.set_number_of_keys(Index(number_of_keys));
    vector_function.set_number_of_cameras(is_uniform_camera_ ? 1 : 2);
    FeatureMapContainer feature_maps(number_of_keys);
    for (size_t i = 0; i < number_of_points; i++)
    {
      feature_maps[i].first = 0;
      feature_maps[i].second = Index(i);
      feature_maps[number_of_points + i].first = 1;
      feature_maps[number_of_points + i].second = Index(i);
    }
    vector_function.set_feature_maps(feature_maps);
    ImageCameraMap image_camera_map;
    image_camera_map.push_back(0);
    image_camera_map.push_back(is_uniform_camera_ ? 0 : 1);
    vector_function.set_image_camera_map(image_camera_map);
    vector_function.intrinsic_computations_mask().set();

    if (fix_mask[FIX_POINTS])
    {
      vector_function.set_fix_points(points_initial);
    }

    Vector3 t_identity = extrinsic_params_identity_initial.rotation() *
                         extrinsic_params_identity_initial.position();
    t_identity *= Scalar(-1);
    Vector3 t_relative = extrinsic_params_relative_initial.rotation() *
                         extrinsic_params_relative_initial.position();
    t_relative *= Scalar(-1);
    BAImageContainer ba_images(2);
    ba_images[0][0] = extrinsic_params_identity_initial.rotation()[0];
    ba_images[0][1] = extrinsic_params_identity_initial.rotation()[1];
    ba_images[0][2] = extrinsic_params_identity_initial.rotation()[2];
    ba_images[0][3] = t_identity[0];
    ba_images[0][4] = t_identity[1];
    ba_images[0][5] = t_identity[2];
    ba_images[1][0] = extrinsic_params_relative_initial.rotation()[0];
    ba_images[1][1] = extrinsic_params_relative_initial.rotation()[1];
    ba_images[1][2] = extrinsic_params_relative_initial.rotation()[2];
    ba_images[1][3] = t_relative[0];
    ba_images[1][4] = t_relative[1];
    ba_images[1][5] = t_relative[2];

    if (fix_mask[FIX_IMAGES])
    {
      vector_function.set_fix_images(ba_images);
    }

    BACameraContainer ba_cameras(is_uniform_camera_ ? 1 : 2);
    Index camera_params_size =
      vector_function.GetIntrinsicParamsSizePerCamera();
    ba_cameras[0].resize(camera_params_size);
    ba_cameras[0][0] = intrinsic_params_identity_initial.k1();
    ba_cameras[0][1] = intrinsic_params_identity_initial.k2();
    ba_cameras[0][2] = intrinsic_params_identity_initial.k3();
    ba_cameras[0][3] = intrinsic_params_identity_initial.d1();
    ba_cameras[0][4] = intrinsic_params_identity_initial.d2();
    ba_cameras[0][5] = intrinsic_params_identity_initial.focal_length();
    ba_cameras[0][6] = intrinsic_params_identity_initial.skew();
    ba_cameras[0][7] = intrinsic_params_identity_initial.principal_point_x();
    ba_cameras[0][8] = intrinsic_params_identity_initial.principal_point_y();
    ba_cameras[0][9] = intrinsic_params_identity_initial.pixel_ratio();
    if (!is_uniform_camera_)
    {
      ba_cameras[1].resize(camera_params_size);
      ba_cameras[1][0] = intrinsic_params_relative_initial.k1();
      ba_cameras[1][1] = intrinsic_params_relative_initial.k2();
      ba_cameras[1][2] = intrinsic_params_relative_initial.k3();
      ba_cameras[1][3] = intrinsic_params_relative_initial.d1();
      ba_cameras[1][4] = intrinsic_params_relative_initial.d2();
      ba_cameras[1][5] = intrinsic_params_relative_initial.focal_length();
      ba_cameras[1][6] = intrinsic_params_relative_initial.skew();
      ba_cameras[1][7] = intrinsic_params_relative_initial.principal_point_x();
      ba_cameras[1][8] = intrinsic_params_relative_initial.principal_point_y();
      ba_cameras[1][9] = intrinsic_params_relative_initial.pixel_ratio();
    }

    if (fix_mask[FIX_CAMERAS])
    {
      vector_function.set_fix_cameras(ba_cameras);
    }
    else
    {
      hs::sfm::ba::CameraConstraint camera_constraint_identity;
      camera_constraint_identity.camera_id = 0;
      camera_constraint_identity.intrinsic_mask.set(
        hs::sfm::ba::INTRINSIC_CONSTRAIN_SKEW);
      camera_constraint_identity.intrinsic_mask.set(
        hs::sfm::ba::INTRINSIC_CONSTRAIN_PIXEL_RATIO);
      vector_function.camera_constraints().push_back(
        camera_constraint_identity);
      if (!is_uniform_camera_)
      {
        hs::sfm::ba::CameraConstraint camera_constraint_relative;
        camera_constraint_relative.camera_id = 1;
        camera_constraint_relative.intrinsic_mask.set(
          hs::sfm::ba::INTRINSIC_CONSTRAIN_SKEW);
        camera_constraint_relative.intrinsic_mask.set(
          hs::sfm::ba::INTRINSIC_CONSTRAIN_PIXEL_RATIO);
        vector_function.camera_constraints().push_back(
          camera_constraint_relative);
      }
    }

    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    XVector x_initial(x_size);
    YVector y_initial(y_size);

    Index x_offset = 0;
    if (!fix_mask[FIX_POINTS])
    {
      auto itr_point = points_initial.begin();
      auto itr_point_end = points_initial.end();
      for (Index i = 0; itr_point != itr_point_end; ++itr_point, i++)
      {
        x_initial.segment(i * VectorFunction::params_per_point_,
                          VectorFunction::params_per_point_) = *itr_point;
      }
      x_offset = vector_function.GetPointParamsSize();
    }

    if (!fix_mask[FIX_IMAGES])
    {
      x_initial.segment(x_offset, VectorFunction::extrinsic_params_per_image_) =
        ba_images[0];
      x_initial.segment(x_offset +
                        VectorFunction::extrinsic_params_per_image_,
                        VectorFunction::extrinsic_params_per_image_) =
        ba_images[1];
      x_offset += vector_function.GetExtrinsicParamsSize();
    }

    if (!fix_mask[FIX_CAMERAS])
    {
      x_initial.segment(x_offset, camera_params_size) = ba_cameras[0];
      if (!is_uniform_camera_)
      {
        x_initial.segment(x_offset + camera_params_size, camera_params_size) =
          ba_cameras[1];
      }
    }

    for (Index i = 0; i < Index(number_of_points); i++)
    {
      y_initial.segment(i * VectorFunction::params_per_key_,
                        VectorFunction::params_per_key_) = keysets[0][i];
      y_initial.segment(number_of_points * VectorFunction::params_per_key_ +
                        i * VectorFunction::params_per_key_,
                        VectorFunction::params_per_key_) = keysets[1][i];
    }

    YCovarianceInverse y_covariance_inverse;
    y_covariance_inverse.SetKeysUniformStdDev(key_stddev_, number_of_keys);
    if (!fix_mask[FIX_CAMERAS])
    {
      Index key_params_size = Index(number_of_points * 4);
      y_initial[key_params_size + 0] = Scalar(0);
      y_initial[key_params_size + 1] = Scalar(1);
      if (!is_uniform_camera_)
      {
        y_initial[key_params_size + 2] = Scalar(0);
        y_initial[key_params_size + 3] = Scalar(1);
      }

      Scalar skew_stddev = Scalar(1e-3);
      Scalar pixel_ratio_stddev = Scalar(1e-2);
      Scalar skew_constraint = Scalar(1) / (skew_stddev * skew_stddev);
      Scalar pixel_ratio_constraint = Scalar(1) / (pixel_ratio_stddev *
                                                  pixel_ratio_stddev);

      y_covariance_inverse.AddConstraint(skew_constraint);
      y_covariance_inverse.AddConstraint(pixel_ratio_constraint);
      if (!is_uniform_camera_)
      {
        y_covariance_inverse.AddConstraint(skew_constraint);
        y_covariance_inverse.AddConstraint(pixel_ratio_constraint);
      }
    }

    XVector x_bundle;
    if (use_ceres)
    {
      CeresOptimizor ceres_optimizor(x_initial);
      if (ceres_optimizor(vector_function, y_initial, y_covariance_inverse,
                          x_bundle) != 0)
      {
        return -1;
      }
    }
    else
    {
      Optimizor optimizor(x_initial, 50,
                          Scalar(1e-3),
                          Scalar(1e-8),
                          Scalar(1e-8));
      if (optimizor(vector_function, y_initial, y_covariance_inverse,
                    x_bundle) != 0)
      {
        return -1;
      }
    }

    x_offset = 0;
    if (!fix_mask[FIX_POINTS])
    {
      points_bundle.resize(number_of_points);
      auto itr_point_bundle = points_bundle.begin();
      auto itr_point_bundle_end = points_bundle.end();
      for (Index i = 0;
           itr_point_bundle != itr_point_bundle_end;
           ++itr_point_bundle, ++i)
      {
        (*itr_point_bundle) =
          x_bundle.segment(i * VectorFunction::params_per_point_,
                           VectorFunction::params_per_point_);
      }
      x_offset = vector_function.GetPointParamsSize();
    }
    else
    {
      points_bundle = points_initial;
    }

    if (!fix_mask[FIX_IMAGES])
    {
      extrinsic_params_identity_bundle.rotation()[0] = x_bundle[x_offset + 0];
      extrinsic_params_identity_bundle.rotation()[1] = x_bundle[x_offset + 1];
      extrinsic_params_identity_bundle.rotation()[2] = x_bundle[x_offset + 2];
      Vector3 c_identity;
      c_identity<<x_bundle[x_offset + 3],
                  x_bundle[x_offset + 4],
                  x_bundle[x_offset + 5];
      extrinsic_params_relative_bundle.rotation()[0] = x_bundle[x_offset + 6];
      extrinsic_params_relative_bundle.rotation()[1] = x_bundle[x_offset + 7];
      extrinsic_params_relative_bundle.rotation()[2] = x_bundle[x_offset + 8];
      Vector3 c_relative;
      c_relative<<x_bundle[x_offset + 9],
                  x_bundle[x_offset + 10],
                  x_bundle[x_offset + 11];

      extrinsic_params_identity_bundle.position() =
        extrinsic_params_identity_bundle.rotation().Inverse() * c_identity;
      extrinsic_params_identity_bundle.position() *= Scalar(-1);
      extrinsic_params_relative_bundle.position() =
        extrinsic_params_relative_bundle.rotation().Inverse() * c_relative;
      extrinsic_params_relative_bundle.position() *= Scalar(-1);

      x_offset += vector_function.GetExtrinsicParamsSize();
    }
    else
    {
      extrinsic_params_identity_bundle = extrinsic_params_identity_initial;
      extrinsic_params_relative_bundle = extrinsic_params_relative_initial;
    }

    if (!fix_mask[FIX_CAMERAS])
    {
      intrinsic_params_identity_bundle.set_k1(
        x_bundle[x_offset + 0]);
      intrinsic_params_identity_bundle.set_k2(
        x_bundle[x_offset + 1]);
      intrinsic_params_identity_bundle.set_k3(
        x_bundle[x_offset + 2]);
      intrinsic_params_identity_bundle.set_d1(
        x_bundle[x_offset + 3]);
      intrinsic_params_identity_bundle.set_d2(
        x_bundle[x_offset + 4]);
      intrinsic_params_identity_bundle.set_focal_length(
        x_bundle[x_offset + 5]);
      intrinsic_params_identity_bundle.set_skew(
        x_bundle[x_offset + 6]);
      intrinsic_params_identity_bundle.set_principal_point_x(
        x_bundle[x_offset + 7]);
      intrinsic_params_identity_bundle.set_principal_point_y(
        x_bundle[x_offset + 8]);
      intrinsic_params_identity_bundle.set_pixel_ratio(
        x_bundle[x_offset + 9]);

      if (is_uniform_camera_)
      {
        intrinsic_params_relative_bundle = intrinsic_params_identity_bundle;
      }
      else
      {
        intrinsic_params_relative_bundle.set_k1(
          x_bundle[x_offset + 10]);
        intrinsic_params_relative_bundle.set_k2(
          x_bundle[x_offset + 11]);
        intrinsic_params_relative_bundle.set_k3(
          x_bundle[x_offset + 12]);
        intrinsic_params_relative_bundle.set_d1(
          x_bundle[x_offset + 13]);
        intrinsic_params_relative_bundle.set_d2(
          x_bundle[x_offset + 14]);
        intrinsic_params_relative_bundle.set_focal_length(
          x_bundle[x_offset + 15]);
        intrinsic_params_relative_bundle.set_skew(
          x_bundle[x_offset + 16]);
        intrinsic_params_relative_bundle.set_principal_point_x(
          x_bundle[x_offset + 17]);
        intrinsic_params_relative_bundle.set_principal_point_y(
          x_bundle[x_offset + 18]);
        intrinsic_params_relative_bundle.set_pixel_ratio(
          x_bundle[x_offset + 19]);
      }
    }
    else
    {
      intrinsic_params_identity_bundle = intrinsic_params_identity_initial;
      intrinsic_params_relative_bundle = intrinsic_params_relative_initial;
    }

    return 0;
  }

private:
  RelativePairGenerator relative_pair_generator_;
  KeysetGenerator keyset_generator_;
  size_t number_of_points_;
  bool is_uniform_camera_;
  Scalar key_stddev_;
};

TEST(TestCameraSharedRelativePair, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestCameraSharedRelativePair<Scalar, ImageDimension>
          Tester;
  typedef Tester::IntrinsicParams IntrinsicParams;

  Scalar pixel_size = 0.0000039;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar focal_length_in_metre_identity = 0.02995452167701055;
  IntrinsicParams intrinsic_params_identity(focal_length_in_metre_identity /
                                            pixel_size,
                                            //0, 0, 0, 1,
                                            0,
                                            -21.669436058,
                                            -44.8644764322,
                                            1,
                                            -0.02529179096221609,
                                            //0, 0, 0, 0
                                            0.23762413973445157,
                                            -0.64208397668697237,
                                            -0.0020605099808780948,
                                            -0.00028706423764766859);
  //IntrinsicParams intrinsic_params_identity(focal_length_in_metre_identity /
  //                                          pixel_size,
  //                                          0, 0, 0, 1, 0, 0, 0, 0, 0);
  Scalar focal_length_in_metre_relative = 0.018858358970276164;
  IntrinsicParams intrinsic_params_relative(focal_length_in_metre_relative /
                                            pixel_size,
                                            //0, 0, 0, 1,
                                            0,
                                            -42.4095312016,
                                            -31.699212823,
                                            1,
                                            -0.0050490462006048831,
                                            //0, 0, 0, 0
                                            0.031293804298609881,
                                            -0.030794820960442223,
                                            -0.00055376548320189127,
                                            -0.00049877717768381476);
  //IntrinsicParams intrinsic_params_relative(focal_length_in_metre_relative /
  //                                          pixel_size,
  //                                          0, 0, 0, 1, 0, 0, 0, 0, 0);
  Scalar ground_resolution = 0.1;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 1.0;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 15;

  size_t number_of_points = 2000;
  bool is_uniform_camera = true;
  Scalar key_stddev = 1.0;

  Tester tester(ground_resolution,
                image_width,
                image_height,
                pixel_size,
                intrinsic_params_identity,
                intrinsic_params_relative,
                lateral_overlap_ratio,
                longitudinal_overlap_ratio,
                scene_max_height,
                camera_height_stddev,
                camera_planar_stddev,
                camera_rotation_stddev,
                number_of_points,
                is_uniform_camera,
                key_stddev);

  ASSERT_EQ(0, tester.Test());
}

}
