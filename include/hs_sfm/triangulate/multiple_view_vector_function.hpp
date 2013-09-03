#ifndef _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_VECTOR_FUNCTION_HPP_
#define _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_VECTOR_FUNCTION_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/utility/cam_type.hpp"

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
  typedef EIGEN_VEC(Scalar, 3) XVec;
  typedef EIGEN_VEC(Scalar, Eigen::Dynamic) YVec;

  typedef hs::sfm::IntrinParam<Scalar> CameraIntrin;
  typedef EIGEN_VECTOR(CameraIntrin) CameraIntrinContainer;
  typedef hs::sfm::ExtrinParam<Scalar> CameraExtrin;
  typedef EIGEN_VECTOR(CameraExtrin) CameraExtrinContainer;

  typedef int Err;

  typedef typename XVec::Index Index;
  static const Index params_per_feature_ = 2;
  static const Index params_per_point_ = 3;

private:
  typedef typename CameraIntrinContainer::size_type SizeType;
  typedef hs::sfm::CamFunc<Scalar> CameraFunction;
  typedef typename CameraFunction::PMat PMatrix;
  typedef EIGEN_VEC(Scalar, params_per_point_ + 1) HomogeneousPoint;
  typedef EIGEN_VEC(Scalar, params_per_feature_ + 1) HomogeneousFeature;

public:
  MultipleViewVectorFunction(){}

  MultipleViewVectorFunction(const CameraIntrinContainer& intrins,
                             const CameraExtrinContainer& extrins)
    : intrins_(intrins), extrins_(extrins) {}

  Err operator()(const XVec& x, YVec& y) const
  {
    Index number_of_camera = Index(intrins_.size());
    if (number_of_camera != Index(extrins_.size()))
    {
      return -1;
    }

    y.resize(number_of_camera * params_per_feature_);

    for (Index i = 0; i < number_of_camera; i++)
    {
      PMatrix p_matrix = CameraFunction::getPMat(intrins_[i], extrins_[i]);

      HomogeneousFeature homogeneous_feature = 
        p_matrix.block(0, 0, 3, 3) * x + 
        p_matrix.block(0, 3, 3, 1);
      homogeneous_feature /= homogeneous_feature[2];

      y.segment(i * params_per_feature_, params_per_feature_) =
        homogeneous_feature.segment(0, params_per_feature_);
    }

    return 0;
  }

  void SetIntrins(const CameraIntrinContainer& intrins)
  {
    intrins_ = intrins;
  }

  void SetExtrins(const CameraExtrinContainer& extrins)
  {
    extrins_ = extrins;
  }

private:
  CameraIntrinContainer intrins_;
  CameraExtrinContainer extrins_;
};

}
}
}

#endif
