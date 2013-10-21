#ifndef _HS_SFM_UTILITY_RANSAC_FIT_PLANE_HPP_
#define _HS_SFM_UTILITY_RANSAC_FIT_PLANE_HPP_

#include "hs_fit/ransac/ransac.hpp"
#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{

template <typename _Point>
struct PointScalarType
{
  typedef typename _Point::Scalar type;
};

//template <typename Scalar>
//struct PtScalarType<Scalar[3]>
//{
//  typedef Scalar type;
//};

//template <typename Scalar>
//struct PtScalarType<Scalar*>
//{
//  typedef Scalar type;
//};

template <typename _Point>
class RansacFitPlane
{
public:
  typedef _Point Point;
  typedef typename hs::fit::PointSetType<Point>::type PointSet;
  typedef typename PointScalarType<Point>::type Scalar;
  typedef EIGEN_VECTOR(Scalar, 4) Plane;
  typedef EIGEN_MATRIX(Scalar, 3, Eigen::Dynamic) Matrix3X;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;

  typedef int Err;

  //计算给定点集的正交回归
  class PlaneOrthogonalRegressionCalculator
  {
  public:
    typedef Plane RansacModel;
    static const size_t model_min_set_size_ = 3;

    int operator()(const PointSet& points, Plane& plane)
    {
      size_t number_of_points = points.size();
      if (number_of_points < model_min_set_size_)
      {
        return -1;
      }

      //pca
      //求均值
      Vector3 mean = Vector3::Zero();
      for (size_t i = 0; i < number_of_points; i++)
      {
        mean[0] += points[i][0];
        mean[1] += points[i][1];
        mean[2] += points[i][2];
      }
      mean /= Scalar(number_of_points);
      Matrix3X A(3, number_of_points);
      for (size_t i = 0; i < number_of_points; i++)
      {
        Vector3 point;
        point << points[i][0],
                 points[i][1],
                 points[i][2];
        A.col(i) = point - mean;
      }
      //协方差矩阵
      Matrix33 covariance = A * A.transpose();
      //求协方差矩阵最小特征值对应的特征向量
      Eigen::EigenSolver<Matrix33> eigen_solver(covariance);
      Scalar min_eigen_value = std::numeric_limits<Scalar>::max();
      int min_id = -1;
      for (int i = 0; i < 3; i++)
      {
        Scalar eigen_value = eigen_solver.eigenvalues()[i].real();
        if (min_eigen_value > eigen_value)
        {
          min_eigen_value = eigen_value;
          min_id = i;
        }
      }
      Vector3 plane_norm = eigen_solver.eigenvectors().col(min_id).real();
      plane_norm /= plane_norm.norm();

      //计算过中值点，且法向量为pca方向的平面
      plane.segment(0, 3) = plane_norm;
      plane[3] = -plane_norm.dot(mean);

      return 0;
    }
  };

  class PointPlaneDistanceCalculator
  {
  public:
    typedef Scalar Distance;

    int operator()(const Point& point, const Plane& plane,
                   Distance& distance)
    {
      Distance norm = plane.segment(0, 3).norm();
      if (norm == 0)
      {
        return -1;
      }

      distance = std::abs(point[0] * plane[0] + 
                          point[1] * plane[1] + 
                          point[2] * plane[2] + plane[3]) / norm;

      return 0;
    }
  };

  typedef hs::fit::Ransac<Point,
                          PlaneOrthogonalRegressionCalculator,
                          PointPlaneDistanceCalculator> FitRansac;
  typedef typename FitRansac::IndexSet IndexSet;

  Err operator()(const PointSet& points, Plane& plane,
                 Scalar distance_threshold) const
  {
    PlaneOrthogonalRegressionCalculator model_calculator;
    PointPlaneDistanceCalculator dist_calculator;
    FitRansac fit_ransac(model_calculator, dist_calculator);

    fit_ransac.SetAlphaThreshold(0.95);
    fit_ransac.SetDistanceThreshold(distance_threshold);
    PointSet refined_points;
    IndexSet inlier_indices;
    PointSet best_points;
    if (fit_ransac(points, refined_points, inlier_indices, best_points) != 0)
    {
      return -1;
    }

    return model_calculator(refined_points, plane);
  }
};

}//sfm
}//hs

#endif
