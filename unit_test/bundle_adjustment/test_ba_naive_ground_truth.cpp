#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/triangulate/multiple_view_levenberg_marquardt_optimizor.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_noised_x_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_noised_y_generator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_synthetic_data_generator.hpp"

namespace
{

template <typename _Scalar>
class TestBANaiveGroundTruth
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef size_t ImageDimension;

private:
  typedef hs::sfm::ba::BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;
  typedef YVector Residuals;
  typedef hs::sfm::ba::BANaiveLevenbergMarquardtOptimizor<VectorFunction>
          Optimizor;
  typedef hs::sfm::ba::BANaiveSyntheticDataGenerator<Scalar, ImageDimension>
          DataGenerator;
  typedef typename DataGenerator::IntrinsicParams IntrinsicParams;
  typedef typename DataGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename DataGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename DataGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename DataGenerator::Point3DContainer Point3DContainer;
  typedef typename DataGenerator::KeysetContainer
          KeysContainer;
  typedef EIGEN_MATRIX(Scalar, 2, 2) Matrix22;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef hs::sfm::ba::BANaiveNoisedXGenerator<Scalar> NoisedXGenerator;
  typedef hs::sfm::ba::BANaiveNoisedYGenerator<Scalar> NoisedYGenerator;
  typedef typename NoisedYGenerator::YCovarianceInverse YCovarianceInverse;
  typedef hs::sfm::triangulate::MultipleViewVectorFunction<Scalar>
          TriangulateVectorFunction;
  typedef typename TriangulateVectorFunction::XVector TriangulateXVector;
  typedef typename TriangulateVectorFunction::YVector TriangulateYVector;
  typedef hs::sfm::triangulate::MultipleViewLevenbergMarquardtOptimizor<
            TriangulateVectorFunction>
          TriangulateOptimizor;
  typedef typename TriangulateOptimizor::YCovarianceInverse
          TriangulateYCovarianceInverse;

public:

  TestBANaiveGroundTruth(Scalar focal_length_in_metre,
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
                         Scalar feature_stddev,
                         Scalar camera_x_rotation_stddev,
                         Scalar camera_y_rotation_stddev,
                         Scalar camera_z_rotation_stddev,
                         Scalar camera_x_pos_stddev,
                         Scalar camera_y_pos_stddev,
                         Scalar camera_z_pos_stddev,
                         Scalar point_x_stddev,
                         Scalar point_y_stddev,
                         Scalar point_z_stddev,
                         size_t number_of_gcp,
                         Scalar gcp_image_stddev)
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
                      north_west_angle),
      feature_stddev_(feature_stddev),
      number_of_gcp_(number_of_gcp),
      gcp_image_stddev_(gcp_image_stddev),
      camera_x_rotation_stddev_(camera_x_rotation_stddev),
      camera_y_rotation_stddev_(camera_y_rotation_stddev),
      camera_z_rotation_stddev_(camera_z_rotation_stddev),
      camera_x_pos_stddev_(camera_x_pos_stddev),
      camera_y_pos_stddev_(camera_y_pos_stddev),
      camera_z_pos_stddev_(camera_z_pos_stddev),
      point_x_stddev_(point_x_stddev),
      point_y_stddev_(point_y_stddev),
      point_z_stddev_(point_z_stddev){}

  Err operator()() const
  {
    VectorFunction vector_function;
    XVector true_x;
    YVector true_y;
    if (data_generator_(vector_function, true_x, true_y) != 0)
    {
      std::cout<<"data generator failed!\n";
      return -1;
    }

    NoisedXGenerator noised_x_generator;
    XVector near_x;
    if (noised_x_generator(vector_function,
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


    YCovarianceInverse y_covariance_inverse;
    GenerateYCovarianceInverse(
      true_y.rows() / VectorFunction::params_per_feature_, y_covariance_inverse);
    NoisedYGenerator noised_y_generator;
    YVector near_y;
    if (noised_y_generator(true_y, y_covariance_inverse, near_y) != 0)
    {
      std::cout<<"noised y vector generator failed!\n";
      return -1;
    }

    Scalar mean_reprojection_error_before =
      CalculateMeanReprojectionError(vector_function, near_x, near_y);
    std::cout<<"mean reprojection error before optimizing:"
             <<mean_reprojection_error_before<<"\n";
    
    Optimizor optimizor(near_x, 100,
                        Scalar(1e-6),
                        Scalar(1e-10),
                        Scalar(1e-10));
    XVector optimized_x;
    if (optimizor(vector_function, near_y, y_covariance_inverse,
                  optimized_x) != 0)
    {
      std::cout<<"optimizor failed!\n";
      return -1;
    }

    Scalar mean_reprojection_error_after =
      CalculateMeanReprojectionError(vector_function, optimized_x, near_y);
    std::cout<<"mean reprojection error after optimizing:"
             <<mean_reprojection_error_after<<"\n";

    //generate gcps
    Point3DContainer gcp_points;
    KeysContainer gcp_image_keys_set;
    hs::sfm::TrackContainer gcp_tracks;
    if (data_generator_.GenerateGCPs(number_of_gcp_, vector_function, true_x,
                                     gcp_points,
                                     gcp_image_keys_set,
                                     gcp_tracks) != 0)
    {
      std::cout<<"generating gcps fails!\n";
      return -1;
    }

    //check accuracy by gcps;
    size_t number_of_gcp_tracks = gcp_tracks.size();
    EIGEN_STD_VECTOR(Vector3) differences;
    for (size_t i = 0; i < number_of_gcp_tracks; i++)
    {
      size_t number_of_views = gcp_tracks[i].size();
      if (number_of_views >= 3)
      {
        IntrinsicParamsContainer intrins(
          number_of_views, data_generator_.GetFocalLengthInPixel());
        ExtrinsicParamsContainer extrins(number_of_views);
        Index triangulate_y_size = 
          number_of_views * TriangulateVectorFunction::params_per_feature_;
        TriangulateYVector true_triangulate_y(triangulate_y_size);
        for (size_t j = 0; j < number_of_views; j++)
        {
          size_t camera_id = gcp_tracks[i][j].first;
          size_t key_id = gcp_tracks[i][j].second;
          extrins[j].rotation()[0] = 
            optimized_x[camera_id * VectorFunction::params_per_camera_ + 0];
          extrins[j].rotation()[1] = 
            optimized_x[camera_id * VectorFunction::params_per_camera_ + 1];
          extrins[j].rotation()[2] = 
            optimized_x[camera_id * VectorFunction::params_per_camera_ + 2];
          Vector3 t = optimized_x.segment(
                        camera_id * VectorFunction::params_per_camera_ + 3, 3);
          extrins[j].position() = -(extrins[j].rotation().Inverse() * t);

          true_triangulate_y.segment(
            j * TriangulateVectorFunction::params_per_feature_,
            TriangulateVectorFunction::params_per_feature_) =
            gcp_image_keys_set[camera_id][key_id];
        }

        TriangulateVectorFunction triangulate_vector_function(
          intrins, extrins);
        TriangulateYCovarianceInverse triangulate_y_covariance =
          TriangulateYCovarianceInverse::Identity(triangulate_y_size,
                                                  triangulate_y_size);
        triangulate_y_covariance *=
          gcp_image_stddev_ * gcp_image_stddev_;
        TriangulateYCovarianceInverse triangulate_y_covariance_inverse =
          TriangulateYCovarianceInverse::Identity(triangulate_y_size,
                                                  triangulate_y_size);
        triangulate_y_covariance_inverse /= 
          gcp_image_stddev_ * gcp_image_stddev_;
        TriangulateYVector noised_triangulate_y;
        hs::math::random::NormalRandomVar<Scalar, Eigen::Dynamic>::Generate(
          true_triangulate_y, triangulate_y_covariance,
          noised_triangulate_y);

        TriangulateOptimizor triangulate_optimizor;
        TriangulateXVector estimated_pt;
        if (triangulate_optimizor(triangulate_vector_function,
                                  noised_triangulate_y,
                                  triangulate_y_covariance_inverse,
                                  estimated_pt) != 0)
        {
          std::cout<<"triangulate optimizor failed!\n";
          return -1;
        }

        Vector3 difference = estimated_pt - gcp_points[i];
        differences.push_back(estimated_pt - gcp_points[i]);
      }
    }

    Vector3 mean = Vector3::Zero();
    Vector3 max = Vector3::Zero();
    for (size_t i = 0; i < differences.size(); i++)
    {
      mean[0] += std::abs(differences[i][0]);
      mean[1] += std::abs(differences[i][1]);
      mean[2] += std::abs(differences[i][2]);

      max[0] = std::max(max[0], std::abs(differences[i][0]));
      max[1] = std::max(max[1], std::abs(differences[i][1]));
      max[2] = std::max(max[2], std::abs(differences[i][2]));
    }
    mean /= Scalar(differences.size());

    if (mean[0] > Scalar(1) ||
        mean[1] > Scalar(1) ||
        mean[2] > Scalar(1))
    {
      return -1;
    }

    return 0;
  }

private:
  void GenerateYCovarianceInverse(
    Index number_of_feature, YCovarianceInverse& y_covariance_inverse) const
  {
    Matrix22 feature_covariance_inverse = Matrix22::Identity();
    Scalar stddev = feature_stddev_ / data_generator_.GetFocalLengthInPixel();
    feature_covariance_inverse /= stddev * stddev;
    y_covariance_inverse.blocks.clear();
    for (Index i = 0; i < number_of_feature; i++)
    {
      y_covariance_inverse.blocks.push_back(feature_covariance_inverse);
    }
  }

  Scalar CalculateMeanReprojectionError(const VectorFunction& vector_function,
                                        const XVector& x,
                                        const YVector& y) const
  {
    YVector noised_y;
    vector_function(x, noised_y);
    if (vector_function(x, noised_y) != 0)
    {
      std::cout<<"vector function failed!\n";
      return Scalar(0);
    }
    Residuals residuals = y - noised_y;

    Index number_of_feature = vector_function.number_of_features();
    Scalar mean_reprojection_error = Scalar(0);
    for (Index i = 0; i < number_of_feature; i++)
    {
      EIGEN_VECTOR(Scalar, 2) residual =
        residuals.segment(i * VectorFunction::params_per_feature_,
                          VectorFunction::params_per_feature_);
      Scalar reprojection_error = residual.norm();
      mean_reprojection_error += reprojection_error;
    }
    mean_reprojection_error /= Scalar(number_of_feature) / 
                               data_generator_.GetFocalLengthInPixel();

    return mean_reprojection_error;
  }

private:
  DataGenerator data_generator_;
  Scalar feature_stddev_;
  size_t number_of_gcp_;
  Scalar gcp_image_stddev_;
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

TEST(TestBANaiveGroundTruth, DoubleTest)
{
  typedef double Scalar;

  typedef TestBANaiveGroundTruth<Scalar> Test;
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

  Scalar camera_x_rotation_stddev = 1;
  Scalar camera_y_rotation_stddev = 1;
  Scalar camera_z_rotation_stddev = 1;
  Scalar camera_x_pos_stddev = 1;
  Scalar camera_y_pos_stddev = 1;
  Scalar camera_z_pos_stddev = 1;
  Scalar point_x_stddev = 1;
  Scalar point_y_stddev = 1;
  Scalar point_z_stddev = 1;

  Scalar feature_stddev = 1.0;

  size_t number_of_gcp = 100;
  Scalar gcp_image_stddev = 2.0;

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
            feature_stddev, 
            camera_x_rotation_stddev,
            camera_y_rotation_stddev,
            camera_z_rotation_stddev,
            camera_x_pos_stddev,
            camera_y_pos_stddev,
            camera_z_pos_stddev,
            point_x_stddev,
            point_y_stddev,
            point_z_stddev,
            number_of_gcp, gcp_image_stddev);

  ASSERT_EQ(0, test());

}

}