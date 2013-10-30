#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_noised_y_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_noised_y_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_noised_x_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_synthetic_data_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_relative_data_generator.hpp"

namespace
{

template <typename _Scalar>
class TestBAGCPConstrainedGauge
{
public:
  public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef size_t ImageDimension;

private:
  typedef hs::sfm::ba::BAGCPConstrainedVectorFunction<Scalar> GCPVectorFunction;
  typedef hs::sfm::ba::BANaiveVectorFunction<Scalar> NaiveVectorFunction;
  typedef typename GCPVectorFunction::XVector XVector;
  typedef typename GCPVectorFunction::YVector YVector;
  typedef typename GCPVectorFunction::Index Index;
  typedef YVector Residuals;
  typedef hs::sfm::ba::BAGCPConstrainedLevenbergMarquardtOptimizor<
            GCPVectorFunction> GCPOptimizor;
  typedef typename GCPOptimizor::YCovarianceInverse
                   GCPYCovarianceInverse;
  typedef hs::sfm::ba::BANaiveLevenbergMarquardtOptimizor<
            NaiveVectorFunction> NaiveOptimizor;
  typedef typename NaiveOptimizor::YCovarianceInverse
                   NaiveYCovarianceInverse;
  typedef hs::sfm::ba::BAGCPConstrainedSyntheticDataGenerator<Scalar,
                                                              ImageDimension>
          DataGenerator;
  typedef typename DataGenerator::IntrinsicParams IntrinsicParams;
  typedef typename DataGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename DataGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename DataGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename DataGenerator::Point3DContainer Point3DContainer;
  typedef typename DataGenerator::KeysetContainer KeysetContainer;
  typedef EIGEN_MATRIX(Scalar, 2, 2) Matrix22;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_MATRIX(Scalar, 4, 4) Matrix44;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef hs::sfm::ba::BANaiveNoisedXGenerator<Scalar> NoisedXGenerator;
  typedef hs::sfm::ba::BAGCPConstrainedNoisedYGenerator<Scalar>
          GCPNoisedYGenerator;
  typedef hs::sfm::ba::BANaiveNoisedYGenerator<Scalar>
          NaiveNoisedYGenerator;
  typedef hs::sfm::ba::BAGCPConstrainedRelativeDataGenerator<Scalar>
          RelativeDataGenerator;
  typedef typename RelativeDataGenerator::Rotation Rotation;
  typedef typename RelativeDataGenerator::Translate Translate;

public:
  TestBAGCPConstrainedGauge(
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
    size_t number_of_gcps,
    Scalar gcp_planar_stddev,
    Scalar gcp_height_stddev,
    Scalar feature_stddev,
    Scalar feature_noise_stddev,
    Scalar camera_x_rotation_stddev,
    Scalar camera_y_rotation_stddev,
    Scalar camera_z_rotation_stddev,
    Scalar camera_x_pos_stddev,
    Scalar camera_y_pos_stddev,
    Scalar camera_z_pos_stddev,
    Scalar point_x_stddev,
    Scalar point_y_stddev,
    Scalar point_z_stddev,
    size_t number_of_check_points,
    Scalar check_point_image_stddev)
    : data_generator_(focal_length_in_metre,
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
                      north_west_angle,
                      number_of_gcps),
      gcp_planar_stddev_(gcp_planar_stddev),
      gcp_height_stddev_(gcp_height_stddev),
      feature_stddev_(feature_stddev),
      feature_noise_stddev_(feature_noise_stddev),
      camera_x_rotation_stddev_(camera_x_rotation_stddev),
      camera_y_rotation_stddev_(camera_y_rotation_stddev),
      camera_z_rotation_stddev_(camera_z_rotation_stddev),
      camera_x_pos_stddev_(camera_x_pos_stddev),
      camera_y_pos_stddev_(camera_y_pos_stddev),
      camera_z_pos_stddev_(camera_z_pos_stddev),
      point_x_stddev_(point_x_stddev),
      point_y_stddev_(point_y_stddev),
      point_z_stddev_(point_z_stddev),
      number_of_check_points_(number_of_check_points),
      check_point_image_stddev_(check_point_image_stddev) {}

  Err operator()() const
  {
    GCPVectorFunction gcp_vector_function;
    XVector true_abs_x;
    YVector gcp_true_abs_y;
    if(data_generator_(gcp_vector_function, true_abs_x, gcp_true_abs_y) != 0)
    {
      std::cout<<"data generator failed!\n";
      return -1;
    }

    NaiveVectorFunction naive_vector_function(gcp_vector_function);
    YVector naive_true_abs_y;
    if (naive_vector_function(true_abs_x, naive_true_abs_y) != 0)
    {
      std::cout<<"naive vector function failed!\n";
      return -1;
    }

    RelativeDataGenerator relative_data_generator;
    Scalar similar_scale;
    Rotation similar_rotation;
    Translate similar_translate;
    XVector true_rel_x;
    if (relative_data_generator(gcp_vector_function, true_abs_x,
                                true_rel_x,
                                similar_scale,
                                similar_rotation,
                                similar_translate) != 0)
    {
      std::cout<<"relative data generator failed!\n";
      return -1;
    }

    YVector naive_true_rel_y = naive_true_abs_y;

    NoisedXGenerator noised_x_generator;
    XVector near_rel_x;
    if (noised_x_generator(gcp_vector_function,
                           true_rel_x,
                           camera_x_rotation_stddev_,
                           camera_y_rotation_stddev_,
                           camera_z_rotation_stddev_,
                           camera_x_pos_stddev_ / similar_scale,
                           camera_y_pos_stddev_ / similar_scale,
                           camera_z_pos_stddev_ / similar_scale,
                           point_x_stddev_ / similar_scale,
                           point_y_stddev_ / similar_scale,
                           point_z_stddev_ / similar_scale,
                           near_rel_x))
    {
      std::cout<<"noised x vector generator failed!\n";
      return -1;
    }

    NaiveYCovarianceInverse naive_y_covariance_inverse;
    GenerateNaiveYCovarianceInverse(
      naive_vector_function.number_of_features(),
      feature_stddev_,
      naive_y_covariance_inverse);

    NaiveYCovarianceInverse naive_y_noise_covariance_inverse;
    GenerateNaiveYCovarianceInverse(
      naive_vector_function.number_of_features(),
      feature_noise_stddev_,
      naive_y_noise_covariance_inverse);
    NaiveNoisedYGenerator naive_noised_y_generator;
    YVector naive_near_rel_y;
    if (naive_noised_y_generator(naive_true_rel_y,
                                 naive_y_noise_covariance_inverse,
                                 naive_near_rel_y) != 0)
    {
      std::cout<<"naive noised y vector generator failed!\n";
      return -1;
    }

    Scalar mean_reprojection_error_before =
      CalculateMeanReprojectionError(naive_vector_function,
                                     near_rel_x,
                                     naive_near_rel_y);
    std::cout<<"mean reprojection error before optimizing:"
             <<mean_reprojection_error_before<<"\n";

    NaiveOptimizor naive_optimizor(near_rel_x, 100,
                                   Scalar(1e-6),
                                   Scalar(1e-10),
                                   Scalar(1e-10));
    XVector naive_optimized_rel_x;
    if (naive_optimizor(naive_vector_function, naive_near_rel_y,
                        naive_y_covariance_inverse,
                        naive_optimized_rel_x) != 0)
    {
      std::cout<<"naive optimizor failed!\n";
      return -1;
    }

    Scalar naive_mean_reprojection_error_after =
      CalculateMeanReprojectionError(naive_vector_function,
                                     naive_optimized_rel_x,
                                     naive_near_rel_y);
    std::cout<<"mean reprojection error after naive optimizing:"
             <<naive_mean_reprojection_error_after<<"\n";

    Scalar estimated_scale;
    Rotation estimated_rotation;
    Translate estimated_translate;
    EstimateSimilarTransform(naive_vector_function,
                             true_abs_x,
                             naive_optimized_rel_x,
                             gcp_vector_function.number_of_gcps(),
                             estimated_scale,
                             estimated_rotation,
                             estimated_translate);

    XVector naive_optimized_abs_x;
    SimilarTransformXVector(naive_vector_function,
                            estimated_scale,
                            estimated_rotation,
                            estimated_translate,
                            naive_optimized_rel_x,
                            naive_optimized_abs_x);
    //SimilarTransformXVector(naive_vector_function,
    //                        similar_scale,
    //                        similar_rotation,
    //                        similar_translate,
    //                        naive_optimized_rel_x,
    //                        naive_optimized_abs_x);

    XVector gcp_near_abs_x = naive_optimized_abs_x;
    GCPYCovarianceInverse gcp_y_covariance_inverse;
    GenerateGCPYCovarianceInverse(
      gcp_vector_function.number_of_features(),
      gcp_vector_function.number_of_gcps(),
      feature_stddev_,
      gcp_y_covariance_inverse);

    GCPYCovarianceInverse gcp_y_noise_covariance_inverse;
    GenerateGCPYCovarianceInverse(
      gcp_vector_function.number_of_features(),
      gcp_vector_function.number_of_gcps(),
      feature_noise_stddev_,
      gcp_y_noise_covariance_inverse);
    GCPNoisedYGenerator gcp_noised_y_generator;
    YVector gcp_near_abs_y;
    if (gcp_noised_y_generator(gcp_true_abs_y,
                               gcp_y_noise_covariance_inverse,
                               gcp_near_abs_y) != 0)
    {
      std::cout<<"gcp noised y vector generator failed!\n";
      return -1;
    }

    Index feature_size = gcp_vector_function.GetYFeatureSize();
    gcp_near_abs_y.segment(0, feature_size) = naive_near_rel_y;

    GCPOptimizor gcp_optimizor(gcp_near_abs_x, 100,
                               Scalar(1e-6),
                               Scalar(1e-10),
                               Scalar(1e-10));

    XVector gcp_optimized_abs_x;
    if (gcp_optimizor(gcp_vector_function, gcp_near_abs_y,
                      gcp_y_covariance_inverse,
                      gcp_optimized_abs_x) != 0)
    {
      std::cout<<"gcp optimizor failed!\n";
      return -1;
    }

    Scalar gcp_mean_reprojection_error_after =
      CalculateMeanReprojectionError(gcp_vector_function,
                                     gcp_optimized_abs_x,
                                     gcp_near_abs_y);
    std::cout<<"mean reprojection error after gcp optimizing:"
             <<gcp_mean_reprojection_error_after<<"\n";

    XVector naive_diff_x = naive_optimized_abs_x - true_abs_x;
    XVector gcp_diff_x = gcp_optimized_abs_x - true_abs_x;

    Index point_id_start = gcp_vector_function.GetCameraParamsSize();
    Index number_of_points = gcp_vector_function.number_of_points();
    Scalar naive_planar_mean = 0;
    Scalar naive_height_mean = 0;
    Scalar naive_planar_max = 0;
    Scalar naive_height_max = 0;
    Scalar planar_threshold = 0.8;
    Scalar height_threshold = 0.8;
    Index naive_number_of_over_planar_threshold = 0;
    Index naive_number_of_over_height_threshold = 0;

    Scalar gcp_planar_mean = 0;
    Scalar gcp_height_mean = 0;
    Scalar gcp_planar_max = 0;
    Scalar gcp_height_max = 0;
    Index gcp_number_of_over_planar_threshold = 0;
    Index gcp_number_of_over_height_threshold = 0;
    for (Index i = 0; i < number_of_points; i++)
    {
      Index point_id =
        point_id_start + i * GCPVectorFunction::params_per_point_;
      Scalar naive_planar_diff = naive_diff_x.segment(point_id, 2).norm();
      Scalar naive_height_diff = std::abs(naive_diff_x[point_id + 2]);
      naive_planar_mean += naive_planar_diff;
      naive_height_mean += naive_height_diff;
      naive_planar_max = std::max(naive_planar_max, naive_planar_diff);
      naive_height_max = std::max(naive_height_max, naive_height_diff);
      if (naive_planar_diff > planar_threshold)
      {
        naive_number_of_over_planar_threshold++;
      }
      if (naive_height_diff > height_threshold)
      {
        naive_number_of_over_height_threshold++;
      }

      Scalar gcp_planar_diff = gcp_diff_x.segment(point_id, 2).norm();
      Scalar gcp_height_diff = std::abs(gcp_diff_x[point_id + 2]);
      gcp_planar_mean += gcp_planar_diff;
      gcp_height_mean += gcp_height_diff;
      gcp_planar_max = std::max(gcp_planar_max, gcp_planar_diff);
      gcp_height_max = std::max(gcp_height_max, gcp_height_diff);
      if (gcp_planar_diff > planar_threshold)
      {
        gcp_number_of_over_planar_threshold++;
      }
      if (gcp_height_diff > height_threshold)
      {
        gcp_number_of_over_height_threshold++;
      }
    }

    naive_planar_mean /= Scalar(number_of_points);
    naive_height_mean /= Scalar(number_of_points);

    gcp_planar_mean /= Scalar(number_of_points);
    gcp_height_mean /= Scalar(number_of_points);

    std::cout<<"naive_planar_mean:"<<naive_planar_mean<<"\n";
    std::cout<<"naive_height_mean:"<<naive_height_mean<<"\n";
    std::cout<<"naive_planar_max:"<<naive_planar_max<<"\n";
    std::cout<<"naive_height_max:"<<naive_height_max<<"\n";
    std::cout<<"naive_number_of_over_planar_threshold:"
             <<naive_number_of_over_planar_threshold<<"\n";
    std::cout<<"naive_number_of_over_height_threshold:"
             <<naive_number_of_over_height_threshold<<"\n";

    std::cout<<"gcp_planar_mean:"<<gcp_planar_mean<<"\n";
    std::cout<<"gcp_height_mean:"<<gcp_height_mean<<"\n";
    std::cout<<"gcp_planar_max:"<<gcp_planar_max<<"\n";
    std::cout<<"gcp_height_max:"<<gcp_height_max<<"\n";
    std::cout<<"gcp_number_of_over_planar_threshold:"
             <<gcp_number_of_over_planar_threshold<<"\n";
    std::cout<<"gcp_number_of_over_height_threshold:"
             <<gcp_number_of_over_height_threshold<<"\n";

    if (gcp_planar_mean < naive_planar_mean &&
        gcp_height_mean < naive_height_mean)
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }

private:
  void GenerateGCPYCovarianceInverse(
    Index number_of_features, Index number_of_gcps, Scalar feature_stddev,
    GCPYCovarianceInverse& gcp_y_covariance_inverse) const
  {
    gcp_y_covariance_inverse.naive_y_covariance_inverse.blocks.clear();
    Index feature_size =
      number_of_features * GCPVectorFunction::params_per_feature_;
    Index y_size = feature_size +
                   number_of_gcps * GCPVectorFunction::params_per_point_;
    Matrix22 feature_covariance_inverse = Matrix22::Identity();
    Scalar stddev = feature_stddev / data_generator_.GetFocalLengthInPixel();
    feature_covariance_inverse /= stddev * stddev;
    Matrix33 gcp_covariance_inverse = Matrix33::Identity();
    gcp_covariance_inverse(0, 0) /= gcp_planar_stddev_ * gcp_planar_stddev_;
    gcp_covariance_inverse(1, 1) /= gcp_planar_stddev_ * gcp_planar_stddev_;
    gcp_covariance_inverse(2, 2) /= gcp_height_stddev_ * gcp_height_stddev_;
    for (Index i = 0; i < number_of_features; i++)
    {
      gcp_y_covariance_inverse.naive_y_covariance_inverse.blocks.push_back(
        feature_covariance_inverse);
    }
    for (Index i = 0; i < number_of_gcps; i++)
    {
      gcp_y_covariance_inverse.gcp_blocks.push_back(gcp_covariance_inverse);
    }
  }

  void GenerateNaiveYCovarianceInverse(
    Index number_of_features, Scalar feature_stddev,
    NaiveYCovarianceInverse& naive_y_covariance_inverse) const
  {
    naive_y_covariance_inverse.blocks.clear();
    Index feature_size =
      number_of_features * NaiveVectorFunction::params_per_feature_;
    Index y_size = feature_size;
    Matrix22 feature_covariance_inverse = Matrix22::Identity();
    Scalar stddev = feature_stddev / data_generator_.GetFocalLengthInPixel();
    feature_covariance_inverse /= stddev * stddev;
    for (Index i = 0; i < number_of_features; i++)
    {
      naive_y_covariance_inverse.blocks.push_back(feature_covariance_inverse);
    }
  }

  template <typename VectorFunction>
  Scalar CalculateMeanReprojectionError(
    const VectorFunction& vector_function,
    const XVector& x,
    const YVector& y) const
  {
    YVector noised_y;
    if (vector_function(x, noised_y) != 0)
    {
      std::cout<<"vector function failed!\n";
      return Scalar(0);
    }

    Residuals residuals = y - noised_y;

    Index number_of_features = vector_function.number_of_features();
    Scalar mean_reprojection_error = Scalar(0);
    for (Index i = 0; i < number_of_features; i++)
    {
      EIGEN_VECTOR(Scalar, 2) residual =
        residuals.segment(i * GCPVectorFunction::params_per_feature_,
                          GCPVectorFunction::params_per_feature_);
      Scalar reprojection_error = residual.norm();
      mean_reprojection_error += reprojection_error;
    }
    mean_reprojection_error /= Scalar(number_of_features) /
                               data_generator_.GetFocalLengthInPixel();

    return mean_reprojection_error;
  }

  Err SimilarTransformXVector(const NaiveVectorFunction& vector_function,
                              Scalar similar_scale,
                              const Rotation& similar_rotation,
                              const Translate& similar_translate,
                              const XVector& rel_x,
                              XVector& abs_x) const
  {
    Index x_size = vector_function.GetXSize();
    abs_x.resize(x_size);
    Index number_of_cameras = vector_function.number_of_cameras();
    for (Index i = 0; i < number_of_cameras; i++)
    {
      Vector3 rel_camera_rotation_vector =
        rel_x.segment(i * NaiveVectorFunction::params_per_camera_, 3);
      Vector3 rel_camera_t =
        rel_x.segment(i * NaiveVectorFunction::params_per_camera_ + 3, 3);
      Rotation rel_camera_rotation(rel_camera_rotation_vector);
      Rotation abs_camera_rotation = rel_camera_rotation *
                                     similar_rotation.Inverse();
      Vector3 abs_camera_t = similar_scale * rel_camera_t;
      abs_camera_t = abs_camera_t -
                     (abs_camera_rotation * similar_translate);

      abs_x[i * NaiveVectorFunction::params_per_camera_ + 0] =
        abs_camera_rotation[0];
      abs_x[i * NaiveVectorFunction::params_per_camera_ + 1] =
        abs_camera_rotation[1];
      abs_x[i * NaiveVectorFunction::params_per_camera_ + 2] =
        abs_camera_rotation[2];
      abs_x.segment(i * NaiveVectorFunction::params_per_camera_ + 3, 3) =
        abs_camera_t;
    }

    Index number_of_points = vector_function.number_of_points();
    Index point_start_id = vector_function.GetCameraParamsSize();
    for (Index i = 0; i < number_of_points; i++)
    {
      Vector3 rel_point =
        rel_x.segment(
          point_start_id + i * NaiveVectorFunction::params_per_point_,
          NaiveVectorFunction::params_per_point_);
      Vector3 abs_point = similar_scale * (similar_rotation * rel_point) +
                          similar_translate;
      abs_x.segment(
        point_start_id + i * NaiveVectorFunction::params_per_point_,
        NaiveVectorFunction::params_per_point_) = abs_point;
    }

    return 0;
  }

  Err EstimateSimilarTransform(const NaiveVectorFunction& vector_function,
                               const XVector& true_abs_x,
                               const XVector& optimized_rel_x,
                               size_t number_of_gcps,
                               Scalar& estimated_scale,
                               Rotation& estimated_rotation,
                               Translate& estimated_translate) const
  {
    EIGEN_STD_VECTOR(Vector3) rel_points(number_of_gcps);
    EIGEN_STD_VECTOR(Vector3) abs_points(number_of_gcps);

    Index gcp_start_id =
      vector_function.GetXSize() -
      number_of_gcps * NaiveVectorFunction::params_per_point_;
    for (size_t i = 0; i < number_of_gcps; i++)
    {
      rel_points[i] =
        optimized_rel_x.segment(gcp_start_id +
        i * NaiveVectorFunction::params_per_point_,
        NaiveVectorFunction::params_per_point_);

      abs_points[i] =
        true_abs_x.segment(gcp_start_id +
                           i * NaiveVectorFunction::params_per_point_,
                           NaiveVectorFunction::params_per_point_);
    }

    EIGEN_STD_VECTOR(Vector3) rel_points_centroid = rel_points;
    EIGEN_STD_VECTOR(Vector3) abs_points_centroid = abs_points;

    Vector3 rel_point_center = Vector3::Zero();
    Vector3 abs_point_center = Vector3::Zero();

    for (size_t i = 0; i < number_of_gcps; i++)
    {
      rel_point_center += rel_points_centroid[i];
      abs_point_center += abs_points_centroid[i];
    }
    rel_point_center /= Scalar(number_of_gcps);
    abs_point_center /= Scalar(number_of_gcps);

    for (size_t i = 0; i < number_of_gcps; i++)
    {
      rel_points_centroid[i] -= rel_point_center;
      abs_points_centroid[i] -= abs_point_center;
    }

    Matrix33 M = Matrix33::Zero();
    for (Index i = 0; i < Index(number_of_gcps); i++)
    {
      M += rel_points_centroid[i] * abs_points_centroid[i].transpose();
    }

    Matrix44 N;
    N << M(0, 0) + M(1, 1) + M(2, 2), M(1, 2) - M(2, 1), M(2, 0) - M(0, 2), M(0, 1) - M(1, 0),
         M(1, 2) - M(2, 1), M(0, 0) - M(1, 1) - M(2, 2), M(0, 1) + M(1, 0), M(2, 0) + M(0, 2),
         M(2, 0) - M(0, 2), M(0, 1) + M(1, 0), M(1, 1) - M(0, 0) - M(2, 2), M(1, 2) + M(2, 1),
         M(0, 1) - M(1, 0), M(2, 0) + M(0, 2), M(1, 2) + M(2, 1), M(2, 2) - M(0, 0) - M(1, 1);

    Eigen::EigenSolver<Matrix44> eigen_solver(N);
    Scalar max_eigen_value = -std::numeric_limits<Scalar>::max();
    int max_id = -1;
    for (Index i = 0; i < 4; i++)
    {
      Scalar eigen_value = eigen_solver.eigenvalues()[i].real();
      if (max_eigen_value < eigen_value)
      {
        max_eigen_value = eigen_value;
        max_id = i;
      }
    }

    Eigen::Quaternion<Scalar, Eigen::AutoAlign> q;
    q.coeffs().segment(0, 3) =
      eigen_solver.eigenvectors().col(max_id).real().segment(1, 3);
    q.w() = eigen_solver.eigenvectors().col(max_id).real()(0);
    Scalar maximizer = 0;
    for (size_t i = 0; i < number_of_gcps; i++)
    {
      Vector3 rel_point_abs = q._transformVector(rel_points_centroid[i]);
      maximizer += rel_point_abs.dot(abs_points_centroid[i]);
    }

    estimated_rotation = q;

    Scalar sum_norm = 0.0;
    for (size_t i = 0; i < number_of_gcps; i++)
    {
      sum_norm += rel_points_centroid[i].squaredNorm();
    }

    estimated_scale = maximizer / sum_norm;

    estimated_translate = abs_point_center -
                          estimated_scale *
                          (estimated_rotation * rel_point_center);
    
    return 0;
  }

private:
  DataGenerator data_generator_;
  Scalar gcp_planar_stddev_;
  Scalar gcp_height_stddev_;
  Scalar feature_stddev_;
  Scalar feature_noise_stddev_;
  size_t number_of_check_points_;
  Scalar check_point_image_stddev_;
  Scalar camera_x_rotation_stddev_;
  Scalar camera_y_rotation_stddev_;
  Scalar camera_z_rotation_stddev_;
  Scalar camera_x_pos_stddev_;
  Scalar camera_y_pos_stddev_;
  Scalar camera_z_pos_stddev_;
  Scalar point_x_stddev_;
  Scalar point_y_stddev_;
  Scalar point_z_stddev_;
};

TEST(TestBAGCPConstrainedGauge, DoubleTest)
{

  typedef double Scalar;

  typedef TestBAGCPConstrainedGauge<Scalar> Test;
  typedef Test::ImageDimension ImageDimension;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 7;
  size_t number_of_cameras_in_strip = 20;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 3000;
  Scalar lateral_overlap_ratio = 0.7;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 100;
  Scalar camera_height_stddev = 5;
  Scalar camera_plannar_stddev = 5;
  Scalar camera_rotation_stddev = 5;
  Scalar north_west_angle = 60; 
  size_t number_of_gcps = 10;

  Scalar gcp_planar_stddev = 0.05;
  Scalar gcp_height_stddev = 0.1;

  Scalar camera_x_rotation_stddev = 1;
  Scalar camera_y_rotation_stddev = 1;
  Scalar camera_z_rotation_stddev = 1;
  Scalar camera_x_pos_stddev = 1;
  Scalar camera_y_pos_stddev = 1;
  Scalar camera_z_pos_stddev = 1;
  Scalar point_x_stddev = 1;
  Scalar point_y_stddev = 1;
  Scalar point_z_stddev = 1;

  Scalar feature_stddev = 1;
  Scalar feature_noise_stddev = 1;

  size_t number_of_check_points = 10;
  Scalar check_point_image_stddev = 2.0;

  Test test(focal_length_in_metre,
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
            north_west_angle,
            number_of_gcps,
            gcp_planar_stddev,
            gcp_height_stddev,
            feature_stddev,
            feature_noise_stddev,
            camera_x_rotation_stddev,
            camera_y_rotation_stddev,
            camera_z_rotation_stddev,
            camera_x_pos_stddev,
            camera_y_pos_stddev,
            camera_z_pos_stddev,
            point_x_stddev,
            point_y_stddev,
            point_z_stddev,
            number_of_check_points,
            check_point_image_stddev);

  ASSERT_EQ(0, test());

}

}