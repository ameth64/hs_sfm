#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_sfm/sfm_utility/similar_transform_estimator.hpp"

namespace
{

template <typename _Scalar>
class TestSimilarTransformEstimator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

private:
  typedef hs::sfm::SimilarTransformEstimator<Scalar> Estimator;
  typedef typename Estimator::Point Point;
  typedef typename Estimator::PointContainer PointContainer;

public:
  typedef typename Estimator::Rotation Rotation;
  typedef typename Estimator::Translate Translate;

public:
  TestSimilarTransformEstimator(
    size_t number_of_points,
    Scalar stddev,
    const Rotation& rotation_prior,
    const Translate& translate_prior,
    Scalar scale_prior)
    : number_of_points_(number_of_points),
      stddev_(stddev),
      rotation_prior_(rotation_prior),
      translate_prior_(translate_prior),
      scale_prior_(scale_prior) {}

public:
  Err Test()
  {
    Point rel_point_min;
    rel_point_min << -1,
                     -1,
                     -1;
    Point rel_point_max = -rel_point_min;

    PointContainer rel_points_true(number_of_points_);
    for (size_t i = 0; i < number_of_points_; i++)
    {
      hs::math::random::UniformRandomVar<Scalar, 3>::Generate(
        rel_point_min, rel_point_max, rel_points_true[i]);
    }

    PointContainer abs_points_true(number_of_points_);
    for (size_t i = 0; i < number_of_points_; i++)
    {
      Point rel_point = rel_points_true[i];
      abs_points_true[i] = scale_prior_ * (rotation_prior_ * rel_point) +
                           translate_prior_;

    }

    typedef EIGEN_MATRIX(Scalar, 3, 3) PointCovariance;
    PointCovariance rel_covariance = PointCovariance::Identity();
    PointCovariance abs_covariance = PointCovariance::Identity();
    rel_covariance *= stddev_ * stddev_ / scale_prior_ / scale_prior_;
    abs_covariance *= stddev_ * stddev_;
    PointContainer rel_points_noised = rel_points_true;
    PointContainer abs_points_noised = abs_points_true;
    for (size_t i = 0; i < number_of_points_; i++)
    {
      hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
        rel_points_true[i], rel_covariance, rel_points_noised[i]);
      hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
        abs_points_true[i], abs_covariance, abs_points_noised[i]);
    }

    Rotation rotation_estimate;
    Translate translate_estimate;
    Scalar scale_estimate;
    Estimator estimator;
    if (estimator(rel_points_noised, abs_points_noised,
                  rotation_estimate,
                  translate_estimate,
                  scale_estimate) != 0)
    {
      return -1;
    }

    PointContainer abs_points_estimate = rel_points_noised;
    for (size_t i = 0; i < number_of_points_; i++)
    {
      Point rel_point = rel_points_noised[i];
      abs_points_estimate[i] =
        scale_estimate * (rotation_estimate * rel_point) + translate_estimate;
        
    }

    Scalar mean_error = Scalar(0);
    for (size_t i = 0; i < number_of_points_; i++)
    {
      Point diff = abs_points_estimate[i] - abs_points_noised[i];
      Scalar error = diff.norm();
      mean_error += error;
    }
    mean_error /= Scalar(number_of_points_);

    Scalar threshold = stddev_ * 2 + 1;
    if (mean_error > threshold)
    {
      return -1;
    }

    return 0;
  }

private:
  size_t number_of_points_;
  Scalar stddev_;
  Rotation rotation_prior_;
  Translate translate_prior_;
  Scalar scale_prior_;
};

TEST(TestSimilarTransformEstimator, SimpleTest)
{
  typedef double Scalar;
  typedef TestSimilarTransformEstimator<Scalar> Tester;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

  typedef Tester::Rotation Rotation;
  typedef Tester::Translate Translate;

  size_t number_of_points = 100;
  Scalar stddev = 0.5;
  Matrix33 rotation_matrix;
  rotation_matrix <<  0.5399779, -0.8415759, 0.0131819,
                      0.8411558,  0.5401282, 0.0268019,
                     -0.0296757, -0.0033843, 0.9995538;
  Rotation rotation_prior(rotation_matrix);
  Translate translate_prior;
  translate_prior << 22.067653670014472,
                     -115.69362180617838,
                     488.11969428744248;
  Scalar scale_prior = 180.83981053672068;

  Tester tester(number_of_points,
                stddev,
                rotation_prior,
                translate_prior,
                scale_prior);
  ASSERT_EQ(0, tester.Test());
}

TEST(TestSimilarTransformEstimator, PriorTest)
{
  typedef double Scalar;
  typedef hs::sfm::SimilarTransformEstimator<Scalar> Estimator;
  typedef Estimator::Point Point;
  typedef Estimator::PointContainer PointContainer;
  typedef Estimator::Rotation Rotation;
  typedef Estimator::Translate Translate;

  PointContainer points_1(4);
  PointContainer points_2(4);

  points_2[0] << 2614304.862,
                 404783.446,
                 116.810;
  points_2[1] << 2611325.079,
                 405029.381,
                 117.931;
  points_2[2] << 2615611.506,
                 411900.269,
                 131.823;
  points_2[3] << 2611160.858,
                 413763.796,
                 108.580;

  points_1[0] << 2614108.822,
                 557851.217,
                 116.810;
  points_1[1] << 2611131.991,
                 558128.369,
                 117.931;
  points_1[2] << 2615490.009,
                 564953.528,
                 131.823;
  points_1[3] << 2611059.331,
                 566863.561,
                 108.580;

  Estimator estimator;
  Rotation rotation_similar;
  Translate translate_similar;
  Scalar scale_similar;
  ASSERT_EQ(0, estimator(points_1, points_2,
                         rotation_similar,
                         translate_similar,
                         scale_similar));

  for (size_t i = 0; i < 4; i++)
  {
    Point point_2_estimate = points_1[i];
    point_2_estimate = rotation_similar * point_2_estimate;
    point_2_estimate = scale_similar * point_2_estimate;
    point_2_estimate = point_2_estimate + translate_similar;
    Point diff = point_2_estimate - points_2[i];
    int bp = 0;
  }

}

}