#ifndef _HS_SFM_SFM_UTILITY_SIMILAR_TRANSFORM_ESTIMATOR_HPP_
#define _HS_SFM_SFM_UTILITY_SIMILAR_TRANSFORM_ESTIMATOR_HPP_

#include <limits>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_math/geometry/rotation.hpp"

namespace hs
{
namespace sfm
{

template <typename _Scalar>
class SimilarTransformEstimator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

public:
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;
  typedef hs::math::geometry::Rotation3D<Scalar> Rotation;
  typedef EIGEN_VECTOR(Scalar, 3) Translate;

private:
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_MATRIX(Scalar, 4, 4) Matrix44;
  typedef typename Point::Index Index;

public:
  Err operator() (const PointContainer& rel_points,
                  const PointContainer& abs_points,
                  Rotation& rotation_similar,
                  Translate& translate_similar,
                  Scalar& scale_similar) const
  {
    size_t number_of_points = rel_points.size();
    if (number_of_points != abs_points.size() ||
        number_of_points < 3)
    {
      return -1;
    }

    Point rel_point_center = Point::Zero();
    Point abs_point_center = Point::Zero();

    for (size_t i = 0; i < number_of_points; i++)
    {
      rel_point_center += rel_points[i];
      abs_point_center += abs_points[i];
    }
    rel_point_center /= Scalar(number_of_points);
    abs_point_center /= Scalar(number_of_points);

    PointContainer rel_points_centroid = rel_points;
    PointContainer abs_points_centroid = abs_points;
    for (size_t i = 0; i < number_of_points; i++)
    {
      rel_points_centroid[i] -= rel_point_center;
      abs_points_centroid[i] -= abs_point_center;
    }

    Matrix33 M = Matrix33::Zero();
    for (Index i = 0; i < Index(number_of_points); i++)
    {
      M += rel_points_centroid[i] * abs_points_centroid[i].transpose();
    }

    Matrix44 N;
    N << M(0,0)+M(1,1)+M(2,2), M(1,2)-M(2,1), M(2,0)-M(0,2), M(0,1)-M(1,0),
         M(1,2)-M(2,1), M(0,0)-M(1,1)-M(2,2), M(0,1)+M(1,0), M(2,0)+M(0,2),
         M(2,0)-M(0,2), M(0,1)+M(1,0), M(1,1)-M(0,0)-M(2,2), M(1,2)+M(2,1),
         M(0,1)-M(1,0), M(2,0)+M(0,2), M(1,2)+M(2,1), M(2,2)-M(0,0)-M(1,1);

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
    for (size_t i = 0; i < number_of_points; i++)
    {
      Point rel_point_abs = q._transformVector(rel_points_centroid[i]);
      maximizer += rel_point_abs.dot(abs_points_centroid[i]);
    }

    rotation_similar = q;

    Scalar sum_norm = 0.0;
    for (size_t i = 0; i < number_of_points; i++)
    {
      sum_norm += rel_points_centroid[i].squaredNorm();
    }

    scale_similar = maximizer / sum_norm;

    translate_similar = abs_point_center -
                        scale_similar *
                        (rotation_similar * rel_point_center);

    return 0;
  }
};

}
}

#endif
