#ifndef _HS_SFM_PROJECTIVE_PMATRIX_DLT_CALCULATOR_HPP_
#define _HS_SFM_PROJECTIVE_PMATRIX_DLT_CALCULATOR_HPP_

#include <utility>

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

template <typename _Scalar>
class PMatrixDLTCalculator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_MATRIX(Scalar, 3, 4) PMatrix;
  typedef std::pair<Key, Point> Correspondence;
  typedef EIGEN_STD_VECTOR(Correspondence) CorrespondenceContainer;

private:
  typedef EIGEN_VECTOR(Scalar, 3) HKey;
  typedef EIGEN_VECTOR(Scalar, 4) HPoint;
  typedef EIGEN_MATRIX(Scalar, 3, 3) KeyTransform;
  typedef EIGEN_MATRIX(Scalar, 4, 4) PointTransform;
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic) ParamMatrix;
  typedef EIGEN_VECTOR(Scalar, 12) PVector;

  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
public:
  Err operator() (const CorrespondenceContainer& correspondences,
                  PMatrix& p_matrix) const
  {
    size_t number_of_correspondences = correspondences.size();
    if (number_of_correspondences < 6)
    {
      return -1;
    }

    KeyTransform key_transform = KeyTransform::Identity();
    PointTransform point_transform = PointTransform::Identity();
    GetNormalTransform(correspondences, key_transform, point_transform);

    ParamMatrix param_matrix(number_of_correspondences * 2, 12);

    for (size_t i = 0; i < number_of_correspondences; i++)
    {
      Key transformed_key = TransformKey(key_transform,
                                         correspondences[i].first);
      Point transformed_point = TransformPoint(point_transform,
                                               correspondences[i].second);
      param_matrix.template block<1, 4>(i * 2 + 0, 0).setZero();
      param_matrix.template block<1, 3>(i * 2 + 0, 4) =
        -transformed_point.transpose();
      param_matrix(i * 2 + 0, 7) = -Scalar(1);
      param_matrix.template block<1, 3>(i * 2 + 0, 8) =
        transformed_key[1] *
        transformed_point.transpose();
      param_matrix(i * 2 + 0, 11) = transformed_key[1];
      param_matrix.template block<1, 3>(i * 2 + 1, 0) =
        transformed_point.transpose();
      param_matrix(i * 2 + 1, 3) = Scalar(1);
      param_matrix.template block<1, 4>(i * 2 + 1, 4).setZero();
      param_matrix.template block<1, 3>(i * 2 + 1, 8) =
        -transformed_key[0] *
        transformed_point.transpose();
      param_matrix(i * 2 + 1, 11) = -transformed_key[0];
    }

    Eigen::JacobiSVD<ParamMatrix> svd(param_matrix, Eigen::ComputeThinV);
    PVector p_vector = svd.matrixV().col(11);

    PMatrix p_matrix_transformed;
    p_matrix_transformed << p_vector[0], p_vector[1], p_vector[2], p_vector[3],
                            p_vector[4], p_vector[5], p_vector[6], p_vector[7],
                            p_vector[8], p_vector[9], p_vector[10],p_vector[11];

    UnNormalizePMatrix(p_matrix_transformed, key_transform, point_transform,
                       p_matrix);

    return 0;
  }

private:
  void GetNormalTransform(const CorrespondenceContainer& correspondences,
                          KeyTransform& key_transform,
                          PointTransform& point_transform) const
  {
    size_t number_of_correspondences = correspondences.size();
    Key key_centroid = Key::Zero();
    Point point_centroid = Point::Zero();
    for (size_t i = 0; i < number_of_correspondences; i++)
    {
      key_centroid += correspondences[i].first;
      point_centroid += correspondences[i].second;
    }
    key_centroid /= Scalar(number_of_correspondences);
    point_centroid /= Scalar(number_of_correspondences);

    Scalar key_mean_distance = Scalar(0);
    Scalar point_mean_distance = Scalar(0);
    for (size_t i = 0; i < number_of_correspondences; i++)
    {
      key_mean_distance +=
        (correspondences[i].first - key_centroid).norm();
      point_mean_distance +=
        (correspondences[i].second - point_centroid).norm();
    }
    key_mean_distance /= Scalar(number_of_correspondences);
    point_mean_distance /= Scalar(number_of_correspondences);
    Scalar key_scale = Scalar(1) / key_mean_distance * std::sqrt(Scalar(2));
    Scalar point_scale = Scalar(1) / point_mean_distance * std::sqrt(Scalar(3));
    key_transform.setIdentity();
    key_transform(0, 0) = key_scale;
    key_transform(1, 1) = key_scale;
    key_transform(0, 2) = -key_scale * key_centroid[0];
    key_transform(1, 2) = -key_scale * key_centroid[1];
    point_transform.setIdentity();
    point_transform.template block<3, 3>(0, 0) *= point_scale;
    point_transform.template block<3, 1>(0, 3) = -point_scale * point_centroid;
  }

  inline Key TransformKey(const KeyTransform& key_transform,
                          const Key& key) const
  {
    return key_transform.template block<2, 2>(0, 0) * key +
           key_transform.template block<2, 1>(0, 2);
  }

  inline Point TransformPoint(const PointTransform& point_transform,
                              const Point& point) const
  {
    return point_transform.template block<3, 3>(0, 0) * point +
           point_transform.template block<3, 1>(0, 3);
  }

  void UnNormalizePMatrix(const PMatrix& p_matrix_transformed,
                          const KeyTransform& key_transform,
                          const PointTransform& point_transform,
                          PMatrix& p_matrix) const
  {
    p_matrix = key_transform.inverse() *
               p_matrix_transformed *
               point_transform;
  }

};

}
}
}

#endif
