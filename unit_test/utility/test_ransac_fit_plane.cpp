#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_sfm/sfm_utility/ransac_fit_plane.hpp"

namespace
{

template <typename _Scalar>
class TestRansacFitPlannar
{
public:
  typedef _Scalar Scalar;
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_VECTOR(Scalar, 4) Vector4;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef Vector3 Point;

  /**
   *  测试ransac拟合三维点平面。
   *  平面真值为xy平面作相似变换后的平面。
   */
  static int Test(const Matrix33& rotation,
                  const Vector3& translate,
                  Scalar scale,
                  const Vector2& min,
                  const Vector2& max,
                  Scalar sigma,
                  Scalar outlier_ratio,
                  Scalar outlier_threshold,
                  size_t number_of_points)
  {
    Matrix33 covariance = Matrix33::Identity();
    covariance *= sigma;
    std::vector<Point> points;
    for (size_t i = 0; i < number_of_points; i++)
    {
      Vector2 planar_point;
      hs::math::random::UniformRandomVar<Scalar, 2>::Generate(
        min, max, planar_point);

      Vector3 mean;
      mean << planar_point[0],
              planar_point[1],
              0;
      Vector3 point;
      hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
        mean, covariance, point);

      Scalar random;
      hs::math::random::UniformRandomVar<Scalar, 1>::Generate(
        Scalar(0), Scalar(1), random);
      if (random < outlier_ratio)
      {
        Scalar outlier;
        hs::math::random::UniformRandomVar<Scalar, 1>::Generate(
          outlier_threshold, outlier_threshold * 3, outlier);

        point[0] += outlier;
        point[1] += outlier;
        point[2] += outlier;
      }

      point = rotation * point * scale + translate;

      Point p;
      p[0] = point[0];
      p[1] = point[1];
      p[2] = point[2];
      points.push_back(p);
    }

    Vector4 plane;
    hs::sfm::RansacFitPlane<Point> ransac_fitter;
    if (ransac_fitter(points, plane, outlier_threshold) != 0) return -1;

    size_t number_of_outliers = 0;
    for (size_t i = 0; i < number_of_points; i++)
    {
      Vector4 point_homogeneous;
      point_homogeneous << points[i][0],
                           points[i][1],
                           points[i][2],
                           1;
      Scalar distance = point_homogeneous.dot(plane) /
                        plane.segment(0, 3).norm();
      if (distance > outlier_threshold * scale)
      {
        number_of_outliers++;
      }
    }

    Scalar statistical_outliers_ratio =
      Scalar(number_of_outliers) / Scalar(number_of_points);

    return statistical_outliers_ratio < outlier_ratio * 2 ? 0 : -1;
  }
};

TEST(TestRansacFitPlannar, SimpleTest)
{
  typedef double Scalar;
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_VECTOR(Scalar, 4) Vector4;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

  Matrix33 rotation = Matrix33::Identity();
  Vector3 translate;
  translate << 0,
               0,
               10;
  Scalar scale = 5.0;

  Vector2 min;
  min << -10,
         -10;
  Vector2 max;
  max << 10,
         10;

  Scalar sigma = 1.0;
  Scalar outlier_ratio = 0.2;
  Scalar outlier_threshold = 5;
  size_t number_of_points = 100;

  ASSERT_EQ(0, TestRansacFitPlannar<Scalar>::
    Test(rotation, 
         translate,
         scale,
         min,
         max,
         sigma,
         outlier_ratio,
         outlier_threshold,
         number_of_points));
}

}
