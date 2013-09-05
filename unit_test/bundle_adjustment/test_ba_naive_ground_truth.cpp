#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/triangulate/multiple_view_levenberg_marquardt_optimizor.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_levenberg_marquardt_optimizor.hpp"

#include "test_ba_naive_noised_x_generator.hpp"
#include "test_ba_naive_noised_y_generator.hpp"
#include "test_ba_naive_base.hpp"

namespace
{

template <typename _Scalar>
class TestBANaiveGroundTruth
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef size_t ImgDim;

private:
  typedef hs::sfm::ba::BANaiveVecFunc<Scalar> VectorFunction;
  typedef typename VectorFunction::XVec XVector;
  typedef typename VectorFunction::YVec YVector;
  typedef typename VectorFunction::Index Index;
  typedef hs::sfm::ba::BANaiveLevenbergMarquardtOptimizor<VectorFunction>
          Optimizor;
  typedef SyntheticDataGenerator<Scalar, ImgDim> DataGenerator;
  typedef typename DataGenerator::Intrin Intrin;
  typedef typename DataGenerator::IntrinContainer IntrinContainer;
  typedef typename DataGenerator::Extrin Extrin;
  typedef typename DataGenerator::ExtrinContainer ExtrinContainer;
  typedef typename DataGenerator::Pt3DContainer Pt3DContainer;
  typedef typename DataGenerator::ImgKeysContainer
          ImageKeysContainer;
  typedef typename DataGenerator::TrackContainer
          TrackContainer;
  typedef EIGEN_MAT(Scalar, 2, 2) Matrix22;
  typedef EIGEN_MAT(Scalar, 3, 3) Matrix33;
  typedef EIGEN_VEC(Scalar, 3) Vector3;
  typedef hs::sfm::ba::BANaiveNoisedXGenerator<Scalar> NoisedXGenerator;
  typedef hs::sfm::ba::BANaiveNoisedYGenerator<Scalar> NoisedYGenerator;
  typedef typename NoisedYGenerator::YCovarianceInverse YCovarianceInverse;
  typedef hs::sfm::triangulate::MultipleViewVectorFunction<Scalar>
          TriangulateVectorFunction;
  typedef typename TriangulateVectorFunction::XVec TriangulateXVector;
  typedef typename TriangulateVectorFunction::YVec TriangulateYVector;
  typedef hs::sfm::triangulate::MultipleViewLevenbergMarquardtOptimizor<
            TriangulateVectorFunction>
          TriangulateOptimizor;
  typedef typename TriangulateOptimizor::YCovarianceInverse
          TriangulateYCovarianceInverse;

public:

  TestBANaiveGroundTruth(Scalar f, size_t strip_num,
                         size_t cams_num_in_strip, Scalar ground_resolution,
                         ImgDim image_width, ImgDim image_height,
                         Scalar pix_size,
                         size_t pts_num,
                         Scalar lateral_overlap,
                         Scalar longitudinal_overlap,
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
    : data_generator_(f, strip_num, cams_num_in_strip, ground_resolution,
                      image_width, image_height, pix_size, pts_num,
                      lateral_overlap, longitudinal_overlap, scene_max_height,
                      camera_height_stddev, camera_plannar_stddev,
                      camera_rotation_stddev, north_west_angle),
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
    GenerateYCovarianceInverse(true_y.rows() / VectorFunction::m_paramsPerFeat,
                               y_covariance_inverse);
    NoisedYGenerator noised_y_generator;
    YVector near_y;
    if (noised_y_generator(true_y, y_covariance_inverse, near_y) != 0)
    {
      std::cout<<"noised y vector generator failed!\n";
      return -1;
    }
    
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

    //generate gcps
    Pt3DContainer gcp_points;
    ImageKeysContainer gcp_image_keys_set;
    TrackContainer gcp_tracks;
    if (data_generator_.genGCPs(number_of_gcp_, vector_function, true_x,
                                gcp_points,
                                gcp_image_keys_set,
                                gcp_tracks) != 0)
    {
      std::cout<<"generating gcps fails!\n";
      return -1;
    }

    //check accuracy by gcps;
    size_t number_of_gcp_tracks = gcp_tracks.size();
    EIGEN_VECTOR(Vector3) differences;
    for (size_t i = 0; i < number_of_gcp_tracks; i++)
    {
      size_t number_of_views = gcp_tracks[i].size();
      if (number_of_views >= 3)
      {
        IntrinContainer intrins(number_of_views, data_generator_.getFInPix());
        ExtrinContainer extrins(number_of_views);
        Index triangulate_y_size = 
          number_of_views * TriangulateVectorFunction::params_per_feature_;
        TriangulateYVector true_triangulate_y(triangulate_y_size);
        for (size_t j = 0; j < number_of_views; j++)
        {
          size_t camera_id = gcp_tracks[i][j].first;
          size_t key_id = gcp_tracks[i][j].second;
          extrins[j].m_r[0] = 
            optimized_x[camera_id * VectorFunction::m_paramsPerCam + 0];
          extrins[j].m_r[1] = 
            optimized_x[camera_id * VectorFunction::m_paramsPerCam + 1];
          extrins[j].m_r[2] = 
            optimized_x[camera_id * VectorFunction::m_paramsPerCam + 2];
          Vector3 t = optimized_x.segment(
                        camera_id * VectorFunction::m_paramsPerCam + 3, 3);
          extrins[j].m_c = -(extrins[j].m_r.inverse() * t);

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
        hs::math::random::NormalRandomVar<Scalar, Eigen::Dynamic>::
          normRandomVar(true_triangulate_y, triangulate_y_covariance,
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
    Scalar stddev = feature_stddev_ / data_generator_.getFInPix();
    feature_covariance_inverse /= stddev * stddev;
    y_covariance_inverse.m_blocks.clear();
    for (Index i = 0; i < number_of_feature; i++)
    {
      y_covariance_inverse.m_blocks.push_back(feature_covariance_inverse);
    }
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
  typedef Test::ImgDim ImgDim;

  Scalar f = 0.019;
  size_t strip_num = 5;
  size_t cams_num_in_strip = 30;
  Scalar ground_resolution = 0.1;
  ImgDim image_width = 6000;
  ImgDim image_height = 4000;
  Scalar pix_size = 0.0000039;
  size_t pts_num = 10000;
  Scalar lateral_overlap = 0.7;
  Scalar longitudinal_overlap = 0.8;
  Scalar scene_max_height = 100;
  Scalar camera_height_stddev = 5;
  Scalar camera_plannar_stddev = 5;
  Scalar camera_rotation_stddev = 5;
  Scalar north_west_angle = 60; 

  Scalar camera_x_rotation_stddev = 0.1;
  Scalar camera_y_rotation_stddev = 0.1;
  Scalar camera_z_rotation_stddev = 0.1;
  Scalar camera_x_pos_stddev = 0.1;
  Scalar camera_y_pos_stddev = 0.1;
  Scalar camera_z_pos_stddev = 0.1;
  Scalar point_x_stddev = 0.1;
  Scalar point_y_stddev = 0.1;
  Scalar point_z_stddev = 0.1;

  Scalar feature_stddev = 1.0;

  size_t number_of_gcp = 100;
  Scalar gcp_image_stddev = 2.0;

  Test test(f, strip_num, cams_num_in_strip, ground_resolution,
            image_width, image_height, pix_size, pts_num,
            lateral_overlap, longitudinal_overlap, scene_max_height,
            camera_height_stddev, camera_plannar_stddev,
            camera_rotation_stddev, north_west_angle,
            camera_x_rotation_stddev,
            camera_y_rotation_stddev,
            camera_z_rotation_stddev,
            camera_x_pos_stddev,
            camera_y_pos_stddev,
            camera_z_pos_stddev,
            point_x_stddev,
            point_y_stddev,
            point_z_stddev,
            feature_stddev, number_of_gcp, gcp_image_stddev);

  ASSERT_EQ(0, test());

}

}