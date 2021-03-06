﻿#ifndef _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_VECTOR_FUNCTION_HPP_
#define _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_VECTOR_FUNCTION_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/projective_functions.hpp"

namespace hs
{
namespace sfm
{
namespace triangulate
{

template <typename _Scalar>
class MultipleViewVectorFunction
{
public:
  typedef _Scalar Scalar;
  typedef EIGEN_VECTOR(Scalar, 3) XVector;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) YVector;

  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;

  typedef int Err;

  typedef typename XVector::Index Index;
  static const Index params_per_feature_ = 2;
  static const Index params_per_point_ = 3;

private:
  typedef hs::sfm::ProjectiveFunctions<Scalar> ProjectiveFunctionsType;
  typedef typename ProjectiveFunctionsType::Key KeyType;

public:
  MultipleViewVectorFunction(){}

  MultipleViewVectorFunction(
    const IntrinsicParamsContainer& intrinsic_params_set,
    const ExtrinsicParamsContainer& extrinsic_params_set)
    : intrinsic_params_set_(intrinsic_params_set),
      extrinsic_params_set_(extrinsic_params_set) {}

  Err operator()(const XVector& x, YVector& y) const
  {
    Index number_of_camera = Index(intrinsic_params_set_.size());
    if (number_of_camera != Index(extrinsic_params_set_.size()))
    {
      return -1;
    }

    y.resize(number_of_camera * params_per_feature_);

    for (Index i = 0; i < number_of_camera; i++)
    {
      KeyType key = ProjectiveFunctionsType::WorldPointProjectToImageKey(
        intrinsic_params_set_[i], extrinsic_params_set_[i], x);

      y.segment(i * params_per_feature_, params_per_feature_) = key;
    }

    return 0;
  }

  void set_intrinsic_params_set(
    const IntrinsicParamsContainer& intrinsic_params_set)
  {
    intrinsic_params_set_ = intrinsic_params_set;
  }

  const IntrinsicParamsContainer& intrinsic_params_set() const
  {
    return intrinsic_params_set_;
  }

  void set_extrinsic_params_set(
    const ExtrinsicParamsContainer& extrinsic_params_set)
  {
    extrinsic_params_set_ = extrinsic_params_set;
  }

  const ExtrinsicParamsContainer& extrinsic_params_set() const
  {
    return extrinsic_params_set_;
  }

private:
  IntrinsicParamsContainer intrinsic_params_set_;
  ExtrinsicParamsContainer extrinsic_params_set_;
};

}
}
}

#endif
