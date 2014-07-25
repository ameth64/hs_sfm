#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/calibrate/planar_calibration_generator.hpp"
#include "hs_sfm/calibrate/planar_calibrator.hpp"

namespace
{

template <typename _Scalar>
class TestPlanarCalibrator
{
public:
  typedef _Scalar Scalar;
  typedef size_t ImageDimension;
  typedef hs::sfm::calibrate::PlanarCalibrationGenerator<Scalar,
                                                         ImageDimension>
          Generator;
  typedef hs::sfm::calibrate::PlanarCalibrator<Scalar> Calibrator;
  typedef typename Calibrator::IntrinsicParams IntrinsicParams;
  typedef typename Calibrator::ExtrinsicParams ExtrinsicParams;
  typedef typename Calibrator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename Calibrator::KeyCovariance KeyCovariance;
  typedef typename Calibrator::PointCovariance PointCovariance;
  typedef typename Calibrator::Key Key;
  typedef typename Calibrator::Point Point;
  typedef typename Calibrator::Correspondence Correspondence;
  typedef typename Calibrator::PatternView PatternView;
  typedef typename Calibrator::PatternViewContainer PatternViewContainer;
private:
  typedef EIGEN_VECTOR(Scalar, 10) IntrinsicVector;

public:
  int Test(const IntrinsicParams& intrinsic_params_true,
           Scalar pattern_grid_size,
           size_t number_of_grid_rows,
           size_t number_of_grid_cols,
           size_t number_of_views,
           ImageDimension image_width,
           ImageDimension image_height,
           Scalar point_x_stddev,
           Scalar point_y_stddev,
           Scalar point_z_stddev,
           Scalar image_x_stddev,
           Scalar image_y_stddev,
           size_t number_of_samples)
  {
    int result = 0;
    while (1)
    {
      Generator generator(intrinsic_params_true,
                          pattern_grid_size,
                          number_of_grid_rows,
                          number_of_grid_cols,
                          number_of_views,
                          image_width,
                          image_height);

      ExtrinsicParamsContainer extrinsic_params_set_true;
      PatternViewContainer pattern_views_true;
      result = generator(extrinsic_params_set_true, pattern_views_true);
      if (result != 0) break;

      KeyCovariance key_covariance = KeyCovariance::Identity();
      key_covariance(0, 0) = image_x_stddev * image_x_stddev;
      key_covariance(1, 1) = image_y_stddev * image_y_stddev;
      PointCovariance point_covariance = PointCovariance::Identity();
      point_covariance(0, 0) *= point_x_stddev;
      point_covariance(1, 1) *= point_y_stddev;
      point_covariance(2, 2) *= point_z_stddev;

      IntrinsicVector intrinsic_vector_mean = IntrinsicVector::Zero();
      IntrinsicVector intrinsic_vector_stddev = IntrinsicVector::Zero();
      IntrinsicVector intrinsic_vector_true =
        VectorizeIntrinsicParams(intrinsic_params_true);
      for (size_t i = 0; i < number_of_samples; i++)
      {
        //随机生成加入噪声的数据
        PatternViewContainer pattern_views_noised = pattern_views_true;
        for (size_t j = 0; j < generator.number_of_views(); j++)
        {
          for (size_t k = 0; k < pattern_views_noised[j].size(); k++)
          {
            hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
              pattern_views_true[j][k].first, key_covariance,
              pattern_views_noised[j][k].first);
            hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
              pattern_views_true[j][k].second, point_covariance,
              pattern_views_noised[j][k].second);

          }// for (size_t k = 0; k < pattern_views_noised[j].size(); k++)
        }// for (size_t j = 0; j < generator.number_of_views(); j++)

        //计算标定后的相机内参数
        Calibrator calibrator;
        IntrinsicParams intrinsic_params_estimate;
        ExtrinsicParamsContainer extrinsic_params_set_estimate;
        result = calibrator(pattern_views_noised,
                            key_covariance,
                            point_covariance,
                            intrinsic_params_estimate,
                            extrinsic_params_set_estimate);
        if (result != 0) return result;

        //统计相机内参数标准差
        IntrinsicVector intrinsic_vector_estimate =
          VectorizeIntrinsicParams(intrinsic_params_estimate);

        IntrinsicVector intrinsic_vector_diff =
          intrinsic_vector_estimate - intrinsic_vector_true;
        intrinsic_vector_mean += intrinsic_vector_diff;
        if (i == 0)
        {
          intrinsic_vector_stddev =
            intrinsic_vector_diff.cwiseProduct(intrinsic_vector_diff);
        }
        else if (i == 1)
        {
          intrinsic_vector_stddev +=
            intrinsic_vector_diff.cwiseProduct(intrinsic_vector_diff);
        }
        else
        {
          intrinsic_vector_stddev *= Scalar(i - 1);
          intrinsic_vector_stddev +=
            intrinsic_vector_diff.cwiseProduct(intrinsic_vector_diff);
          intrinsic_vector_stddev /= Scalar(i);
        }
      }// for (size_t i = 0; i < number_of_samples; i++)

      intrinsic_vector_stddev = intrinsic_vector_stddev.cwiseSqrt();

      if (intrinsic_vector_stddev[0] > Scalar(1) ||
          intrinsic_vector_stddev[1] > Scalar(1e-3) ||
          intrinsic_vector_stddev[2] > Scalar(1) ||
          intrinsic_vector_stddev[3] > Scalar(1) ||
          intrinsic_vector_stddev[4] > Scalar(1e-1) ||
          intrinsic_vector_stddev[5] > Scalar(1e-3) ||
          intrinsic_vector_stddev[6] > Scalar(1e-3) ||
          intrinsic_vector_stddev[7] > Scalar(1e-3) ||
          intrinsic_vector_stddev[8] > Scalar(1e-3) ||
          intrinsic_vector_stddev[9] > Scalar(1e-3))
      {
        std::cout<<"intrincic params standard deviation:\n"
                 <<"focal length: "<<intrinsic_vector_stddev[0]<<"\n"
                 <<"skew: "<<intrinsic_vector_stddev[1]<<"\n"
                 <<"principal point x: "<<intrinsic_vector_stddev[2]<<"\n"
                 <<"principal point y: "<<intrinsic_vector_stddev[3]<<"\n"
                 <<"pixel ratio: "<<intrinsic_vector_stddev[4]<<"\n"
                 <<"k1: "<<intrinsic_vector_stddev[5]<<"\n"
                 <<"k2: "<<intrinsic_vector_stddev[6]<<"\n"
                 <<"k3: "<<intrinsic_vector_stddev[7]<<"\n"
                 <<"d1: "<<intrinsic_vector_stddev[8]<<"\n"
                 <<"d2: "<<intrinsic_vector_stddev[9]<<"\n";
      }

      break;
    }
    return result;
  }

private:
  static IntrinsicVector VectorizeIntrinsicParams(
    const IntrinsicParams& intrinsic_params)
  {
    IntrinsicVector intrinsic_vector;
    intrinsic_vector[0] = intrinsic_params.focal_length();
    intrinsic_vector[1] = intrinsic_params.skew();
    intrinsic_vector[2] = intrinsic_params.principal_point_x();
    intrinsic_vector[3] = intrinsic_params.principal_point_y();
    intrinsic_vector[4] = intrinsic_params.pixel_ratio();
    intrinsic_vector[5] = intrinsic_params.k1();
    intrinsic_vector[6] = intrinsic_params.k2();
    intrinsic_vector[7] = intrinsic_params.k3();
    intrinsic_vector[8] = intrinsic_params.d1();
    intrinsic_vector[9] = intrinsic_params.d2();

    return intrinsic_vector;
  }
};

TEST(TestPlanarCalibrator, SimpleTest)
{
  typedef double Scalar;
  typedef TestPlanarCalibrator<Scalar> Tester;
  typedef Tester::IntrinsicParams IntrinsicParams;

  IntrinsicParams intrinsic_params_true(Scalar(3303.44),
                                        Scalar(0),
                                        Scalar(2433.06),
                                        Scalar(1602.12),
                                        Scalar(1),
                                        Scalar(-0.06608472),
                                        Scalar(-0.08491661),
                                        Scalar( 0.00987291),
                                        Scalar(-0.00112938),
                                        Scalar(-0.00013884));
  Scalar pattern_grid_size = Scalar(0.03);
  size_t number_of_grid_rows = 9;
  size_t number_of_grid_cols = 12;
  size_t number_of_views = 10;
  size_t image_width = 4912;
  size_t image_height = 3264;
  Scalar point_x_stddev = Scalar(1e-4);
  Scalar point_y_stddev = Scalar(1e-4);
  Scalar point_z_stddev = Scalar(2e-4);
  Scalar image_x_stddev = Scalar(0.5);
  Scalar image_y_stddev = Scalar(0.5);
  size_t number_of_samples = 1000;

  Tester tester;
  ASSERT_EQ(0, tester.Test(intrinsic_params_true,
                           pattern_grid_size,
                           number_of_grid_rows,
                           number_of_grid_cols,
                           number_of_views,
                           image_width,
                           image_height,
                           point_x_stddev,
                           point_y_stddev,
                           point_z_stddev,
                           image_x_stddev,
                           image_y_stddev,
                           number_of_samples));
}

}
