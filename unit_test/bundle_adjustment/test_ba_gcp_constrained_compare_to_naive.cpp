#include <iostream>
#include <fstream>

#include <gtest/gtest.h>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_noised_y_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_noised_y_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_noised_x_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_synthetic_data_generator.hpp"

namespace
{

template <typename _Scalar>
class TestBAGCPConstrainedCompareToNaive
{
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
  typedef typename DataGenerator::KeysContainer KeysContainer;
  typedef typename DataGenerator::TrackContainer TrackContainer;
  typedef EIGEN_MATRIX(Scalar, 2, 2) Matrix22;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef hs::sfm::ba::BANaiveNoisedXGenerator<Scalar> NoisedXGenerator;
  typedef hs::sfm::ba::BAGCPConstrainedNoisedYGenerator<Scalar>
          GCPNoisedYGenerator;
  typedef hs::sfm::ba::BANaiveNoisedYGenerator<Scalar>
          NaiveNoisedYGenerator;

public:
  TestBAGCPConstrainedCompareToNaive(
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
    XVector true_x;
    YVector gcp_true_y;
    if (data_generator_(gcp_vector_function, true_x, gcp_true_y) != 0)
    {
      std::cout<<"data generator failed!\n";
      return -1;
    }
    NaiveVectorFunction naive_vector_function(gcp_vector_function);
    YVector naive_true_y;
    if (naive_vector_function(true_x, naive_true_y) != 0)
    {
      std::cout<<"naive vector function failed!\n";
      return -1;
    }

    NoisedXGenerator noised_x_generator;
    XVector near_x;
    if (noised_x_generator(gcp_vector_function,
                           true_x,
                           camera_x_rotation_stddev_,
                           camera_y_rotation_stddev_,
                           camera_z_rotation_stddev_,
                           camera_x_pos_stddev_,
                           camera_y_pos_stddev_,
                           camera_z_pos_stddev_,
                           point_x_stddev_,
                           point_y_stddev_,
                           point_z_stddev_,
                           near_x) != 0)
    {
      std::cout<<"noised x vector generator failed!\n";
      return -1;
    }

    GCPYCovarianceInverse gcp_y_covariance_inverse;
    GenerateGCPYCovarianceInverse(
      gcp_vector_function.number_of_features(),
      gcp_vector_function.number_of_gcps(),
      feature_stddev_,
      gcp_y_covariance_inverse);

    NaiveYCovarianceInverse naive_y_covariance_inverse;
    GenerateNaiveYCovarianceInverse(
      naive_vector_function.number_of_features(),
      feature_stddev_,
      naive_y_covariance_inverse);

    GCPYCovarianceInverse gcp_y_noise_covariance_inverse;
    GenerateGCPYCovarianceInverse(
      gcp_vector_function.number_of_features(),
      gcp_vector_function.number_of_gcps(),
      feature_noise_stddev_,
      gcp_y_noise_covariance_inverse);
    GCPNoisedYGenerator gcp_noised_y_generator;
    YVector gcp_near_y;
    if (gcp_noised_y_generator(gcp_true_y,
                               gcp_y_noise_covariance_inverse,
                               gcp_near_y) != 0)
    {
      std::cout<<"gcp noised y vector generator failed!\n";
      return -1;
    }

    NaiveYCovarianceInverse naive_y_noise_covariance_inverse;
    GenerateNaiveYCovarianceInverse(
      naive_vector_function.number_of_features(),
      feature_noise_stddev_,
      naive_y_noise_covariance_inverse);
    NaiveNoisedYGenerator naive_noised_y_generator;
    YVector naive_near_y;
    if (naive_noised_y_generator(naive_true_y,
                                 naive_y_noise_covariance_inverse,
                                 naive_near_y) != 0)
    {
      std::cout<<"naive noised y vector generator failed!\n";
      return -1;
    }

    Index feature_size = gcp_vector_function.GetYFeatureSize();
    //gcp_near_y.segment(0, feature_size) = naive_near_y;
    naive_near_y = gcp_near_y.segment(0, feature_size);

    Scalar gcp_mean_reprojection_error_before =
      CalculateMeanReprojectionError(gcp_vector_function, near_x, gcp_near_y);
    std::cout<<"gcp mean reprojection error before optimizing:"
             <<gcp_mean_reprojection_error_before<<"\n";
    
    GCPOptimizor gcp_optimizor(near_x, 100,
                        Scalar(1e-6),
                        Scalar(1e-10),
                        Scalar(1e-10));
    XVector gcp_optimized_x;
    if (gcp_optimizor(gcp_vector_function, gcp_near_y, gcp_y_covariance_inverse,
                  gcp_optimized_x) != 0)
    {
      std::cout<<"gcp optimizor failed!\n";
      return -1;
    }

    Scalar mean_reprojection_error_after =
      CalculateMeanReprojectionError(gcp_vector_function,
                                     gcp_optimized_x, gcp_near_y);
    std::cout<<"mean reprojection error after gcp optimizing:"
             <<mean_reprojection_error_after<<"\n";

    XVector gcp_diff_x = gcp_optimized_x - true_x;
    
    //naive optimization
    NaiveOptimizor naive_optimizor(near_x, 100,
                        Scalar(1e-6),
                        Scalar(1e-10),
                        Scalar(1e-10));
    XVector naive_optimized_x;
    if (naive_optimizor(naive_vector_function, naive_near_y,
                        naive_y_covariance_inverse,
                        naive_optimized_x) != 0)
    {
      std::cout<<"naive optimizor failed!\n";
      return -1;
    }

    Scalar naive_mean_reprojection_error_after =
      CalculateMeanReprojectionError(naive_vector_function,
                                     naive_optimized_x,
                                     naive_near_y);
    std::cout<<"mean reprojection error after naive optimizing:"
             <<naive_mean_reprojection_error_after<<"\n";

    XVector naive_diff_x = naive_optimized_x - true_x;

    std::ofstream tmp_file("tmp.txt", std::ios::out);
    for (Index i = 0; i < naive_diff_x.rows(); i++)
    {
      tmp_file<<i<<" "<<std::abs(gcp_diff_x[i])
                 <<" "<<std::abs(naive_diff_x[i])
                 <<" "<<std::abs(gcp_diff_x[i]) -
                        std::abs(naive_diff_x[i])<<"\n";
    }

    Scalar gcp_planar_mean = 0;
    Scalar gcp_height_mean = 0;
    Scalar naive_planar_mean = 0;
    Scalar naive_height_mean = 0;
    Index point_id_start = gcp_vector_function.GetCameraParamsSize();
    Index number_of_points = gcp_vector_function.number_of_points();
    for (Index i = 0; i < number_of_points; i++)
    {
      Index point_id =
        point_id_start + i * GCPVectorFunction::params_per_point_;
      gcp_planar_mean += gcp_diff_x.segment(point_id, 2).norm();
      gcp_height_mean += std::abs(gcp_diff_x[point_id + 2]);
      naive_planar_mean += naive_diff_x.segment(point_id, 2).norm();
      naive_height_mean += std::abs(naive_diff_x[point_id + 2]);
    }

    gcp_planar_mean /= Scalar(number_of_points);
    gcp_height_mean /= Scalar(number_of_points);
    naive_planar_mean /= Scalar(number_of_points);
    naive_height_mean /= Scalar(number_of_points);

    tmp_file<<"gcp_planar_mean:"<<gcp_planar_mean<<"\n";
    tmp_file<<"gcp_height_mean:"<<gcp_height_mean<<"\n";
    tmp_file<<"naive_planar_mean:"<<naive_planar_mean<<"\n";
    tmp_file<<"naive_height_mean:"<<naive_height_mean<<"\n";

    return 0;
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

TEST(TestBAGCPConstrainedCompareToNaive, DoubleTest)
{
  typedef double Scalar;

  typedef TestBAGCPConstrainedCompareToNaive<Scalar> Test;
  typedef Test::ImageDimension ImageDimension;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 5;
  size_t number_of_cameras_in_strip = 20;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 2000;
  Scalar lateral_overlap_ratio = 0.7;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 100;
  Scalar camera_height_stddev = 5;
  Scalar camera_plannar_stddev = 5;
  Scalar camera_rotation_stddev = 5;
  Scalar north_west_angle = 60; 
  size_t number_of_gcps = 5;

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