#ifndef _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_VECTOR_FUNCION_HPP_
#define _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_VECTOR_FUNCION_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/sfm_utility/projective_functions.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"

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
  typedef hs::sfm::ProjectiveFunctions<Scalar> ProjectiveFunctionsType;
  typedef typename ProjectiveFunctionsType::IntrinsicParams IntrinsicParams;
  typedef typename ProjectiveFunctionsType::ExtrinsicParams ExtrinsicParams;
  typedef typename ProjectiveFunctionsType::Key Key;

public:
  Err operator() (const XVector& x, YVector& y) const
  {
    Index y_size = GetYSize();
    y.resize(y_size);
    Index number_of_points = Index(points_.size());
    IntrinsicParams intrinsic_params(x[6], x[7], x[8], x[9], x[10],
                                     x[11], x[12], x[13], x[14], x[15]);
    ExtrinsicParams extrinsic_params;
    extrinsic_params.rotation()[0] = x[0];
    extrinsic_params.rotation()[1] = x[1];
    extrinsic_params.rotation()[2] = x[2];
    extrinsic_params.position()[0] = x[3];
    extrinsic_params.position()[1] = x[4];
    extrinsic_params.position()[2] = x[5];
    for (Index i = 0; i < number_of_points; i++)
    {
      const Vector3& point = points_[i];
      Key key = ProjectiveFunctionsType::WorldPointProjectToImageKey(
                  intrinsic_params, extrinsic_params, point);
      y.template segment<2>(i * 2) = key;
    }

    y.template segment<10>(2 * number_of_points) = x.template segment<10>(6);

    return 0;
  }

  Index GetXSize() const
  {
    return 16;
  }

  Index GetYSize() const
  {
    return Index(points_.size() * 2 + 10);
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
