#ifndef _HS_SFM_FUNDAMENTAL_LINEAR_8_POINTS_CALCULATOR_HPP_
#define _HS_SFM_FUNDAMENTAL_LINEAR_8_POINTS_CALCULATOR_HPP_

#include <cmath>

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace fundamental
{

template <typename _Scalar>
class Linear8PointsCalculator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef EIGEN_STD_VECTOR(Key) KeyContainer;
  typedef EIGEN_MATRIX(Scalar, 3, 3) FMatrix;

private:
  typedef EIGEN_MATRIX(Scalar, 3, 3) AffineMatrix; 
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic) ConstraintMatrix;
  typedef typename Key::Index Index;
  typedef EIGEN_VECTOR(Scalar, 9) FVector;

public:
  /**
   *  给定左影像点集以及匹配的右影像点集，使用八点法计算基本矩阵矩阵\f$\mathbf{F}\f$。
   *  
   *  对于左影像点\f$\mathbf{k_l}\f$以及匹配的右影像点\f$\mathbf{k_r}\f$，
   *  计算的基本矩阵\f$\mathbf{F}\f$满足：
   *  \f[
   *    \mathbf{k_r}^T\mathbf{F}\mathbf{k_l} = 0
   *  \f]
   */
  Err operator() (const KeyContainer& keys_left,
                  const KeyContainer& keys_right,
                  FMatrix& f_matrix) const
  {
    size_t number_of_keys = keys_left.size();
    if (number_of_keys != keys_right.size() ||
        number_of_keys < 8)
    {
      return -1;
    }

    AffineMatrix transform_left;
    GetKeysNormalTransform(keys_left, transform_left);
    AffineMatrix transform_right;
    GetKeysNormalTransform(keys_right, transform_right);

    //构造由各对匹配点组成的约束矩阵
    ConstraintMatrix C(Index(number_of_keys), 9);
    for (Index i = 0; i < Index(number_of_keys); i++)
    {
      Key left_normalized_key = TransformKey(transform_left, keys_left[i]);
      Key right_normalized_key = TransformKey(transform_right, keys_right[i]);
      C(i, 0) = right_normalized_key[0] * left_normalized_key[0];
      C(i, 1) = right_normalized_key[0] * left_normalized_key[1];
      C(i, 2) = right_normalized_key[0];
      C(i, 3) = right_normalized_key[1] * left_normalized_key[0];
      C(i, 4) = right_normalized_key[1] * left_normalized_key[1];
      C(i, 5) = right_normalized_key[1];
      C(i, 6) = left_normalized_key[0];
      C(i, 7) = left_normalized_key[1];
      C(i, 8) = Scalar(1);
    }

    Eigen::JacobiSVD<ConstraintMatrix> svd(C, Eigen::ComputeThinV);
    FVector f_vector = svd.matrixV().col(8);
    FMatrix f_matrix_transformed_full_rank;
    f_matrix_transformed_full_rank(0, 0) = f_vector[0];
    f_matrix_transformed_full_rank(0, 1) = f_vector[1];
    f_matrix_transformed_full_rank(0, 2) = f_vector[2];
    f_matrix_transformed_full_rank(1, 0) = f_vector[3];
    f_matrix_transformed_full_rank(1, 1) = f_vector[4];
    f_matrix_transformed_full_rank(1, 2) = f_vector[5];
    f_matrix_transformed_full_rank(2, 0) = f_vector[6];
    f_matrix_transformed_full_rank(2, 1) = f_vector[7];
    f_matrix_transformed_full_rank(2, 2) = f_vector[8];

    FMatrix f_matrix_transformed_rank2;
    NearestRank2Matrix(f_matrix_transformed_full_rank,
                       f_matrix_transformed_rank2);

    UnNormalizeFMatrix(transform_left, transform_right,
                       f_matrix_transformed_rank2,
                       f_matrix);

    return 0;
  }

private:
  void GetKeysNormalTransform(const KeyContainer& keys,
                     AffineMatrix& transform) const
  {
    size_t number_of_keys = keys.size();
    Key centroid = Key::Zero();
    for (size_t i = 0; i < number_of_keys; i++)
    {
      centroid += keys[i];
    }
    centroid /= Scalar(number_of_keys);

    Scalar mean_distance = Scalar(0);
    for (size_t i = 0; i < number_of_keys; i++)
    {
      mean_distance += (keys[i] - centroid).norm();
    }
    mean_distance /= Scalar(number_of_keys);
    Scalar scale = Scalar(1) / mean_distance * std::sqrt(Scalar(2.0));
    transform.setIdentity();
    transform(0, 0) = scale;
    transform(1, 1) = scale;
    transform(0, 2) = -scale * centroid[0];
    transform(1, 2) = -scale * centroid[1];
  }

  inline Key TransformKey(const AffineMatrix& transform,
                          const Key& key) const
  {
    return transform.template block<2, 2>(0, 0) * key +
           transform.template block<2, 1>(0, 2);
  }

  inline void UnNormalizeFMatrix(const AffineMatrix& transform_left,
                          const AffineMatrix& transform_right,
                          const FMatrix& f_matrix_transformed,
                          FMatrix& f_matrix) const
  {
    f_matrix = transform_right.transpose() *
               f_matrix_transformed *
               transform_left;
  }

  void NearestRank2Matrix(const FMatrix full_rank_matrix,
                          FMatrix& rank_2_matrix) const
  {
    Eigen::JacobiSVD<FMatrix> svd(full_rank_matrix,
                                  Eigen::ComputeFullU |
                                  Eigen::ComputeFullV);
    FMatrix diagonal = FMatrix::Zero();
    diagonal(0, 0) = svd.singularValues()[0];
    diagonal(1, 1) = svd.singularValues()[1];
    rank_2_matrix = svd.matrixU() * diagonal * svd.matrixV().transpose();
  }
};

}
}
}

#endif