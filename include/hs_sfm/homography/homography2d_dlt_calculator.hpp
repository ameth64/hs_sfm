#ifndef _HS_SFM_HOMOGRAPHY_HOMOGRAPHY2D_DLT_CALCULATPR_HPP_
#define _HS_SFM_HOMOGRAPHY_HOMOGRAPHY2D_DLT_CALCULATPR_HPP_

#include <utility>

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace homography
{

template <typename _Scalar>
class Homography2DDLTCalculator
{
public:
  typedef _Scalar Scalar;

  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef EIGEN_MATRIX(Scalar, 3, 3) HMatrix;
  typedef std::pair<Key, Key> Correspondence;
  typedef EIGEN_STD_VECTOR(Correspondence) CorrespondenceContainer;

private:
  typedef EIGEN_VECTOR(Scalar, 3) HKey;
  typedef EIGEN_VECTOR(Scalar, 9) HVector;
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic) ParamMatrix;
  typedef EIGEN_MATRIX(Scalar, 3, 3) KeyTransform;

public:
  int operator() (const CorrespondenceContainer& correspondences,
                  HMatrix& h_matrix) const
  {
    size_t number_of_correspondences = correspondences.size();
    if (number_of_correspondences < 4)
    {
      return -1;
    }
    KeyTransform key1_transform = KeyTransform::Identity();
    KeyTransform key2_transform = KeyTransform::Identity();
    //GetNormalTransform(correspondences, key1_transform, key2_transform);

    ParamMatrix param_matrix(number_of_correspondences * 2, 9);
    for (size_t i = 0; i < number_of_correspondences; i++)
    {
      Key transformed_key1 = TransformKey(key1_transform,
                                          correspondences[i].first);
      Key transformed_key2 = TransformKey(key2_transform,
                                          correspondences[i].second);
      param_matrix.template block<1, 3>(i * 2 + 0, 0).setZero();
      param_matrix.template block<1, 2>(i * 2 + 0, 3) =
        -transformed_key1.transpose();
      param_matrix(i * 2 + 0, 5) = Scalar(-1);
      param_matrix.template block<1, 2>(i * 2 + 0, 6) =
        transformed_key2[1] * transformed_key1.transpose();
      param_matrix(i * 2 + 0, 8) = transformed_key2[1];
      param_matrix.template block<1, 2>(i * 2 + 1, 0) =
        transformed_key1.transpose();
      param_matrix(i * 2 + 1, 2) = Scalar(1);
      param_matrix.template block<1, 3>(i * 2 + 1, 3).setZero();
      param_matrix.template block<1, 2>(i * 2 + 1, 6) =
        -transformed_key2[0] * transformed_key1.transpose();
      param_matrix(i * 2 + 1, 8) = -transformed_key2[0];
    }
    Eigen::JacobiSVD<ParamMatrix> svd(param_matrix, Eigen::ComputeFullV);
    HVector h_vector = svd.matrixV().col(svd.matrixV().cols() - 1);

    HMatrix h_matrix_transformed;
    h_matrix_transformed << h_vector[0], h_vector[1], h_vector[2],
                            h_vector[3], h_vector[4], h_vector[5],
                            h_vector[6], h_vector[7], h_vector[8];
    UnnormalizeHMatrix(h_matrix_transformed, key1_transform, key2_transform,
                       h_matrix);
    std::cout<<"param_matrix:\n"<<param_matrix<<"\n";
    std::cout<<"h_matrix:\n"<<h_matrix<<"\n";
    return 0;
  }
private:
  void GetNormalTransform(const CorrespondenceContainer& correspondences,
                          KeyTransform& key1_transform,
                          KeyTransform& key2_transform) const
  {
    size_t number_of_correspondences = correspondences.size();
    Key key1_centroid = Key::Zero();
    Key key2_centroid = Key::Zero();
    for (size_t i = 0; i < number_of_correspondences; i++)
    {
      key1_centroid += correspondences[i].first;
      key2_centroid += correspondences[i].second;
    }
    key1_centroid /= Scalar(number_of_correspondences);
    key2_centroid /= Scalar(number_of_correspondences);

    Scalar key1_mean_distance = Scalar(0);
    Scalar key2_mean_distance = Scalar(0);
    for (size_t i = 0; i < number_of_correspondences; i++)
    {
      key1_mean_distance += (correspondences[i].first - key1_centroid).norm();
      key2_mean_distance += (correspondences[i].second - key2_centroid).norm();
    }
    key1_mean_distance /= Scalar(number_of_correspondences);
    key2_mean_distance /= Scalar(number_of_correspondences);
    Scalar key1_scale = Scalar(1) / key1_mean_distance * std::sqrt(Scalar(2));
    Scalar key2_scale = Scalar(1) / key2_mean_distance * std::sqrt(Scalar(2));
    key1_transform.setIdentity();
    key1_transform(0, 0) = key1_scale;
    key1_transform(1, 1) = key1_scale;
    key1_transform(0, 2) = -key1_scale * key1_centroid[0];
    key1_transform(1, 2) = -key1_scale * key1_centroid[1];
    key2_transform.setIdentity();
    key2_transform(0, 0) = key2_scale;
    key2_transform(1, 1) = key2_scale;
    key2_transform(0, 2) = -key2_scale * key2_centroid[0];
    key2_transform(1, 2) = -key2_scale * key2_centroid[1];
  }

  inline Key TransformKey(const KeyTransform& key_transform,
                          const Key& key) const
  {
    return key_transform.template block<2, 2>(0, 0) * key +
           key_transform.template block<2, 1>(0, 2);
  }

  void UnnormalizeHMatrix(const HMatrix& h_matrix_transformed,
                          const KeyTransform& key1_transform,
                          const KeyTransform& key2_transform,
                          HMatrix& h_matrix) const
  {
    h_matrix = key2_transform.inverse() *
               h_matrix_transformed *
               key1_transform;
  }
};

}
}
}

#endif
