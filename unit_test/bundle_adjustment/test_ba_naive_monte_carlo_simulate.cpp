#include <iostream>

#include <gtest/gtest.h>

#include "hs_test_utility/test_monte_carlo/normal_mle_simulator.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"

#include "test_ba_naive_base.hpp"
#include "test_ba_naive_normal_mle_meta.hpp"

namespace
{

template <typename _Scalar>
class TestBANaiveMonteCarloSimulate
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::ba::BANaiveVecFunc<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef hs::test::mc::NormalMLESimulator<VectorFunction> Simulator;
  typedef typename Simulator::NoisedYVectorGenerator
    NoisedYVectorGenerator;
  typedef typename Simulator::AnalyticXCovarianceCalculator
    AnalyticXCovarianceCalculator;
  typedef typename Simulator::StatisticalXCovarianceCalculator
    StatisticalXCovarianceCalculator;
  typedef typename Simulator::XVectorOptimizor
    XVectorOptimizor;
  typedef typename Simulator::XVector XVector;
  typedef typename Simulator::YVector YVector;
  typedef typename Simulator::YCovarianceInverse YCovarianceInverse;
  typedef typename Simulator::XCovariance XCovariance;

  typedef EIGEN_MAT(Scalar, 2, 2) Mat22;

  static Err Test(const VectorFunction& vector_function,
                  const XVector& true_x,
                  const Mat22& feature_covariance,
                  Scalar x_rotation_stddev,
                  Scalar y_rotation_stddev,
                  Scalar z_rotation_stddev,
                  Scalar x_pos_stddev,
                  Scalar y_pos_stddev,
                  Scalar z_pos_stddev,
                  size_t simulator_size)
  {
    YVector true_y;
    if (vector_function(true_x, true_y) != 0)
    {
      std::cout<<"vector function failed!\n";
      return -1;
    }

    XVector near_x = true_x;
    //put some noise to near_x
    XVectorOptimizor x_optimizor(near_x, 100,
                                 Scalar(1e-6),
                                 Scalar(1e-10),
                                 Scalar(1e-10));

    NoisedYVectorGenerator noised_y_generator;

    AnalyticXCovarianceCalculator analytic_x_covariance_calculator;

    Index number_of_feature = vector_function.getFeatNum();
    YCovarianceInverse y_covariance_inverse;
    if (GenFeatCovInv(feature_covariance, number_of_feature,
                      y_covariance_inverse) != 0)
    {
      std::cout<<"GenFeatCovInv Failed.\n";
      return -1;
    }

    StatisticalXCovarianceCalculator statistical_x_covariance_calculator;
    XCovariance analytic_x_covariance;
    XCovariance statistical_x_covariance;

    Simulator simulator;
    if (simulator(vector_function,
                  x_optimizor,
                  noised_y_generator,
                  analytic_x_covariance_calculator,
                  true_y,
                  y_covariance_inverse,
                  true_x,
                  simulator_size,
                  statistical_x_covariance_calculator,
                  analytic_x_covariance,
                  statistical_x_covariance) != 0)
    {
      std::cout<<"simulator failed!\n";
      return -1;
    }

    const Scalar threshold = Scalar(1e-2);
    Err result = 0;
    for (Index i = 0; i < analytic_x_covariance.rows(); i++)
    {
      for (Index j = 0; j < analytic_x_covariance.cols(); j++)
      {
        Scalar analytic_value = analytic_x_covariance(i, j);
        Scalar statistical_value = statistical_x_covariance(i, j);
        Scalar error = std::abs(analytic_value - statistical_value);
        if (error > threshold)
        {
          std::cout<<"difference between analytic value statistical value "
                     "is to large!\n";
          std::cout<<"analytic_x_covariance["<<i<<", "<<j<<"]="
                   <<analytic_value<<" \n";
          std::cout<<"statistical_x_covariance["<<i<<", "<<j<<"]="
                   <<statistical_value<<" \n";
          result = -1;
        }
      }
    }

    return result;
  }

private:
  static Err GenFeatCovInv(const Mat22& feat_cov,
                           Index number_of_feature,
                           YCovarianceInverse& y_cov_inv)
  {
    y_cov_inv.m_blocks.clear();
    for (Index i = 0; i < number_of_feature; i++)
    {
      y_cov_inv.m_blocks.push_back(feat_cov.inverse());
    }

    return 0;
  }

};

TEST(TestBANaiveMonteCarloSimulate, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImgDim;

  typedef SyntheticDataGenerator<Scalar, ImgDim> DataGen;
  typedef TestBANaiveMonteCarloSimulate<Scalar> Test;
  typedef Test::Mat22 Mat22;
  typedef Test::VectorFunction VectorFunction;
  typedef Test::XVector XVector;
  typedef Test::YVector YVector;

  Scalar f = 0.006;
  size_t strip_num = 1;
  size_t cams_num_in_strip = 4;
  Scalar grd_res = 0.1;
  ImgDim img_w = 6000;
  ImgDim img_h = 4000;
  Scalar pix_size = 0.00000203311408298266;
  size_t pts_num = 10;
  Scalar lateral_overlap = 0.99;
  Scalar longitudinal_overlap = 0.99;
  Scalar scene_max_height = 50;
  Scalar cam_height_dev = 0.01;
  Scalar cam_plannar_dev = 0.01;
  Scalar cam_rot_dev = 0.01;
  Scalar nw_angle = 60;

  VectorFunction ba_vec_func;
  XVector x;
  YVector y;

  DataGen data_gen(f, strip_num, cams_num_in_strip, grd_res,
    img_w, img_h, pix_size, pts_num,
    lateral_overlap, longitudinal_overlap,
    scene_max_height,
    cam_height_dev, cam_plannar_dev, cam_rot_dev, nw_angle);
  ASSERT_EQ(0, data_gen(ba_vec_func, x, y));
  Scalar f_in_pix = f / pix_size;

  Mat22 feat_cov = Mat22::Identity();
  feat_cov *= Scalar(1) / (f_in_pix * f_in_pix);

  ASSERT_EQ(0, Test::Test(ba_vec_func, x, feat_cov,
                          10.0, 10.0, 10.0, 5.0, 5.0, 5.0,
                          100));
}

}