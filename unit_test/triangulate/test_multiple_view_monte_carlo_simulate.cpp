#include <iostream>

#include <gtest/gtest.h>

#include "hs_test_utility/test_monte_carlo/normal_mle_simulator.hpp"
#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/triangulate/multiple_view_vector_function.hpp"

#include "test_multiple_view_base.hpp"
#include "test_multiple_view_normal_mle_meta.hpp"

namespace
{

template <typename _Scalar>
class TestMultipleViewMonteCarloSimulate
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::triangulate::MultipleViewVectorFunction<Scalar>
          VectorFunction;
  typedef hs::test::mc::NormalMLESimulator<VectorFunction> Simulator;
  typedef typename Simulator::NoisedYVectorGenerator
    NoisedYVectorGenerator;
  typedef typename Simulator::AnalyticXCovarianceCalculator
    AnalyticXCovarianceCalculator;
  typedef typename Simulator::StatisticalXCovarianceCalculator
    StatisticalXCovarianceCalculator;
  typedef typename Simulator::XVectorOptimizor
    XVectorOptimizor;
  typedef typename Simulator::ResidualsCalculator
    ResidualsCalculator;
  typedef typename Simulator::MahalanobisDistanceCalculator
    MahalanobisDistanceCalculator;
  typedef typename Simulator::XVector XVector;
  typedef typename Simulator::YVector YVector;
  typedef typename Simulator::YCovarianceInverse YCovarianceInverse;
  typedef typename Simulator::XCovariance XCovariance;

private:
  typedef hs::math::random::NormalRandomVar<Scalar,
                                            VectorFunction::params_per_point_>
          NoisedXVectorGenerator;
  typedef typename XVector::Index Index;
public:

  static Err Test(const VectorFunction& vector_function,
                  const XVector& true_x,
                  Scalar feature_stddev,
                  Scalar point_stddev,
                  size_t simulator_size)
  {
    YVector true_y;
    if (vector_function(true_x, true_y) != 0)
    {
      std::cout<<"vector function failed!\n";
      return -1;
    }

    //generate noised x
    XVector near_x = true_x;
    XCovariance x_covariance = XCovariance::Identity();
    x_covariance *= point_stddev * point_stddev;
    if (NoisedXVectorGenerator::Generate(true_x, x_covariance, near_x))
    {
      std::cout<<"normal random variant generation failed!\n";
      return -1;
    }

    //Monte Carlo simulate
    XVectorOptimizor x_optimizor;
    NoisedYVectorGenerator noised_y_generator;
    AnalyticXCovarianceCalculator analytic_x_covariance_calculator;
    ResidualsCalculator residuals_calculator;
    MahalanobisDistanceCalculator mahalanobis_distance_calculator;
    YCovarianceInverse y_covariance_inverse =
      YCovarianceInverse::Identity(true_y.rows(), true_y.rows());
    y_covariance_inverse /= feature_stddev * feature_stddev;
    StatisticalXCovarianceCalculator statistical_x_covariance_calculator;
    XCovariance analytic_x_covariance;
    XCovariance statistical_x_covariance;
    Scalar residual_mean;

    Simulator simulator;
    if (simulator(vector_function,
                  x_optimizor,
                  residuals_calculator,
                  mahalanobis_distance_calculator,
                  noised_y_generator,
                  analytic_x_covariance_calculator,
                  true_y,
                  y_covariance_inverse,
                  true_x,
                  simulator_size,
                  statistical_x_covariance_calculator,
                  residual_mean,
                  analytic_x_covariance,
                  statistical_x_covariance) != 0)
    {
      std::cout<<"simulator failed!\n";
      return -1;
    }
    residual_mean = std::sqrt(residual_mean / Scalar(true_y.rows()));

    Scalar expected_residual_mean =
      std::sqrt(Scalar(1) - Scalar(3) / Scalar(true_y.rows()));

    if (std::abs(residual_mean - expected_residual_mean) > Scalar(5e-2))
    {
      std::cout<<"Difference between residual mean and expected residual mean "
                 "too large!\n";
      std::cout<<"residual_mean = "<<residual_mean<<".\n";
      std::cout<<"expected_residual_mean = "<<expected_residual_mean<<".\n";
      return -1;
    }

    Err result = 0;
    const Scalar precision = Scalar(1e-1);
    for (Index i = 0; i < VectorFunction::params_per_point_; i++)
    {
      for (Index j = 0; j < VectorFunction::params_per_point_; j++)
      {
        Scalar analytic_value = analytic_x_covariance(i, j);
        Scalar statistical_value = statistical_x_covariance(i, j);
        if (std::abs(analytic_value - statistical_value) > precision)
        {
          std::cout<<"Difference between analytic value and statistical value "
                     "too large!\n";
          std::cout<<"analytic_x_covariance["<<i<<", "<<j<<"] = "
                   <<analytic_value<<".\n";
          std::cout<<"statistical_x_covariance["<<i<<", "<<j<<"] = "
                   <<statistical_value<<".\n";

          result = -1;
        }
      }
    }

    return result;
  }
};

TEST(TestMultipleViewMonteCarloSimulate, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImgDim;
  typedef hs::sfm::triangulate::MultipleViewSyntheicDataGenerator<Scalar,
                                                                  ImgDim>
          DataGenerator;
  typedef DataGenerator::IntrinsicParams IntrinsicParams;
  typedef DataGenerator::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef TestMultipleViewMonteCarloSimulate<Scalar> Test;
  typedef Test::VectorFunction VectorFunction;
  typedef VectorFunction::XVector XVector;
  typedef VectorFunction::YVector YVector;

  Scalar focal_length_in_metre = 0.035;
  size_t number_of_strips = 2;
  size_t number_of_cameras_in_strip = 3;
  Scalar ground_resolution = 0.1;
  ImgDim image_width = 6000;
  ImgDim image_height = 4000;
  Scalar pixel_size = 0.00000444313;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 1;
  Scalar north_west_angle = 60;

  IntrinsicParamsContainer intrinsic_params_set;
  IntrinsicParams intrinsic_params0(8506.46,
                                    0,
                                    2704.86,
                                    1762.62,
                                    1,
                                    0.0599996,
                                    -0.23131,
                                    -0.00460478,
                                    0.000162892,
                                    -0.000169231);
  intrinsic_params_set.push_back(intrinsic_params0);
  intrinsic_params_set.push_back(intrinsic_params0);
  intrinsic_params_set.push_back(intrinsic_params0);
  IntrinsicParams intrinsic_params1(8501.47,
                                    0,
                                    2752.73,
                                    1862.18,
                                    1,
                                    0.0576605,
                                    -0.254137,
                                    0.14565,
                                    -0.00114618,
                                    -0.000143283);
  intrinsic_params_set.push_back(intrinsic_params1);
  intrinsic_params_set.push_back(intrinsic_params1);
  intrinsic_params_set.push_back(intrinsic_params1);

  VectorFunction vector_function;
  XVector x;
  YVector y;

  DataGenerator data_generator(focal_length_in_metre,
                               number_of_strips,
                               number_of_cameras_in_strip,
                               ground_resolution,
                               image_width,
                               image_height,
                               pixel_size,
                               lateral_overlap_ratio,
                               longitudinal_overlap_ratio,
                               scene_max_height,
                               camera_height_stddev,
                               camera_planar_stddev,
                               camera_rotation_stddev,
                               north_west_angle,
                               intrinsic_params_set);
  ASSERT_EQ(0, data_generator(vector_function, x, y));
  ASSERT_EQ(0, Test::Test(vector_function, x, 2.0, 2.0, 100));
}
}
