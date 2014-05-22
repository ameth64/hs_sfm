#ifndef _HS_SFM_ESSENTIAL_EMATRIX_EXTRINSIC_PARAMS_POINTS_CALCULATOR_HPP_
#define _HS_SFM_ESSENTIAL_EMATRIX_EXTRINSIC_PARAMS_POINTS_CALCULATOR_HPP_

#include <utility>

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/triangulate/multiple_view_vector_function.hpp"
#include "hs_sfm/triangulate/multiple_view_levenberg_marquardt_optimizor.hpp"

namespace hs
{
namespace sfm
{
namespace essential
{

template <typename _Scalar>
class EMatrixExtrinsicParamsPointsCalculator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;
  typedef EIGEN_VECTOR(Scalar, 3) HKey;
  typedef std::pair<HKey, HKey> HKeyPair;
  typedef EIGEN_STD_VECTOR(HKeyPair) HKeyPairContainer;
  typedef EIGEN_MATRIX(Scalar, 3, 3) EMatrix;

private:
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;

public:
  Err operator() (const EMatrix& e_matrix,
                  const HKeyPairContainer& key_pairs,
                  ExtrinsicParams& extrin_params,
                  PointContainer& points) const
  {
    Err result = -1;
    EIGEN_MATRIX(Scalar, 3, 3) W;
    W << 0, -1, 0,
         1, 0, 0,
         0, 0, 1;
    EIGEN_MATRIX(Scalar, 3, 3) WT = W.transpose();

    Eigen::JacobiSVD<EMatrix> svd(e_matrix, Eigen::ComputeFullU |
      Eigen::ComputeFullV);
    RMatrix Ra = svd.matrixU() * W * (svd.matrixV().transpose());
    RMatrix Rb = svd.matrixU() * WT * (svd.matrixV().transpose());
    Point t = svd.matrixU().col(2);
    if (Ra.determinant() < 0) Ra *= -1;
    if (Rb.determinant() < 0) Rb *= -1;

    if (TriangulatePoints(key_pairs, Rb, t, points) != 0)
    {
      return -1;
    }
    if (IsValidPoints(Rb, t, points) == 0)
    {
      extrin_params.rotation() = Rb;
      extrin_params.position() = -Rb.transpose() * t;
      return 0;
    }

    if (TriangulatePoints(key_pairs, Rb, -t, points) != 0)
    {
      return -1;
    }
    if (IsValidPoints(Rb, -t, points) == 0)
    {
      extrin_params.rotation() = Rb;
      extrin_params.position() = Rb.transpose() * t;
      return 0;
    }

    if (TriangulatePoints(key_pairs, Ra, t, points) != 0)
    {
      return -1;
    }
    if (IsValidPoints(Ra, t, points) == 0)
    {
      extrin_params.rotation() = Ra;
      extrin_params.position() = -Ra.transpose() * t;
      return 0;
    }

    if (TriangulatePoints(key_pairs, Ra, -t, points) != 0)
    {
      return -1;
    }
    if (IsValidPoints(Ra, -t, points) == 0)
    {
      extrin_params.rotation() = Ra;
      extrin_params.position() = Ra.transpose() * t;
      return 0;
    }

    return -1;
  }

private:
  Err TriangulatePoints(const HKeyPairContainer& key_pairs,
                        const RMatrix& R,
                        const Point& t,
                        PointContainer& points) const
  {
    typedef hs::sfm::triangulate::MultipleViewVectorFunction<Scalar>
            VectorFunction;
    typedef typename VectorFunction::XVector XVector;
    typedef typename VectorFunction::YVector YVector;
    typedef typename VectorFunction::Index Index;
    typedef hs::sfm::triangulate::MultipleViewLevenbergMarquardtOptimizor<
              VectorFunction> Optimizor;
    typedef typename Optimizor::YCovarianceInverse YCovarianceInverse;
    typedef typename ExtrinsicParams::Rotation Rotation;
    typedef typename ExtrinsicParams::Position Position;
    IntrinsicParamsContainer intrinsic_params_pair;
    ExtrinsicParamsContainer extrinsic_params_pair;
    IntrinsicParams identity_intirnsic_params(1, 0, 0, 0, 1);
    intrinsic_params_pair.push_back(identity_intirnsic_params);
    intrinsic_params_pair.push_back(identity_intirnsic_params);

    ExtrinsicParams identity_extrinsic_params;
    extrinsic_params_pair.push_back(identity_extrinsic_params);
    Position position = -R.transpose() * t;
    Rotation rotation = R;
    ExtrinsicParams relative_extrinsic_params(rotation, position);
    extrinsic_params_pair.push_back(relative_extrinsic_params);
    VectorFunction vector_function(intrinsic_params_pair,
                                   extrinsic_params_pair);

    size_t number_of_key_pairs = key_pairs.size();
    YCovarianceInverse y_covariance_inverse =
      YCovarianceInverse::Identity(4, 4);
    points.resize(number_of_key_pairs);
    for (size_t i = 0; i < number_of_key_pairs; i++)
    {
      const HKey& key_left = key_pairs[i].first;
      const HKey& key_right = key_pairs[i].second;
      YVector y(4);
      y[0] = key_left[0] / key_left[2];
      y[1] = key_left[1] / key_left[2];
      y[2] = key_right[0] / key_right[2];
      y[3] = key_right[1] / key_right[2];

      XVector x;
      Optimizor optimizor;
      if (optimizor(vector_function, y, y_covariance_inverse, x) != 0)
      {
        return -1;
      }

      points[i] = x;
    }

    return 0;
  }

  Err IsValidPoints(const RMatrix& R,
                    const Point& t,
                    const PointContainer& points) const
  {
    typedef EIGEN_MATRIX(Scalar, 3, 4) PMatrix;
    PMatrix P_left;
    P_left << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0;
    PMatrix P_right;
    P_right.block(0, 0, 3, 3) = R;
    P_right.block(0, 3, 3, 1) = t;

    size_t front_left = 0;
    size_t behind_left = 0;
    size_t front_right = 0;
    size_t behind_right = 0;
    for (size_t i = 0; i < points.size(); i++)
    {
      Scalar z_left = (P_left.block(2, 0, 1, 3) * points[i])[0] +
                      P_left(2, 3);
      Scalar z_right = (P_right.block(2, 0, 1, 3) * points[i])[0] +
                       P_right(2, 3);
      if (z_left < 0) front_left++; else behind_left++;
      if (z_right < 0) front_right++; else behind_right++;
    }

    if (front_left > behind_left && front_right > behind_right)
    {
      return 0;
    }
    else
    {
      return -1;
    }

    return -1;
  }
};

}
}
}

#endif
