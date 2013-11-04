#ifndef _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_VECTOR_FUNCION_HPP_
#define _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_VECTOR_FUNCION_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

template <typename _Scalar>
class IntrinsicConstrainedVectorFunction
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) XVector;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) YVector;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;
  typedef typename XVector::Index Index;

private:
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;

public:
  Err operator() (const XVector& x, YVector& y) const
  {
    Index y_size = GetYSize();
    y.resize(y_size);
    Index number_of_points = Index(points_.size());
    Vector3 rotation  = x.template segment<3>(0);
    Vector3 translate =  x.template segment<3>(3);
    Matrix33 K;
    K << x[6], x[7], x[8],
         0, x[6] * x[10], x[9],
         0, 0, 1;
    Scalar theta = rotation.norm();
    if (theta != Scalar(0))
    {
      rotation /= theta;
    }
    for (Index i = 0; i < number_of_points; i++)
    {
      const Vector3& point = points_[i];
      Vector3 key_homogeneous;
      if (theta == Scalar(0))
      {
        key_homogeneous = K * (point + translate);
      }
      else
      {
        key_homogeneous =
          cos(theta) * point +
          sin(theta) * rotation.cross(point) +
          (1 - cos(theta)) * point.dot(rotation) * rotation +
          translate;
        key_homogeneous = K * key_homogeneous;
      }
      key_homogeneous /= key_homogeneous[2];
      y.template segment<2>(i * 2) = key_homogeneous.template segment<2>(0);
    }

    y.template segment<5>(2 * number_of_points) = x.template segment<5>(6);

    return 0;
  }

  Index GetXSize() const
  {
    return 11;
  }

  Index GetYSize() const
  {
    return Index(points_.size() * 2 + 5);
  }

  const Point3DContainer&  points() const
  {
    return points_;
  }

  Point3DContainer& points()
  {
    return points_;
  }

private:
  Point3DContainer points_;
};

}
}
}

#endif
