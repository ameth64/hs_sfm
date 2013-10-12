#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_JACOBIAN_MATRIX_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_JACOBIAN_MATRIX_HPP_

#include "hs_sfm/bundle_adjustment/ba_naive_jacobian_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar,
          typename _Index,
          _Index params_per_feature,
          _Index params_per_camera,
          _Index params_per_point>
class BAGCPConstrainedJacobianMatrix
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef _Index Index;
  typedef BANaiveJacobianMatrix<Scalar, Index, params_per_feature,
                                               params_per_camera,
                                               params_per_point>
          NaiveJacobianMatrix;

public:
  inline NaiveJacobianMatrix& naive_jacobian_matrix()
  {
    return naive_jacobian_matrix_;
  }

  inline const NaiveJacobianMatrix& naive_jacobian_matrix() const
  {
    return naive_jacobian_matrix_;
  }

  inline Index number_of_gcps() const
  {
    return number_of_gcps_;
  }

  inline void set_number_of_gcps(Index number_of_gcps)
  {
    number_of_gcps_ = number_of_gcps;
  }

  Scalar coeff(Index i, Index j) const
  {
    Index y_feature_size = 
      naive_jacobian_matrix_.camera_derivatives().size() * params_per_feature;
    Index naive_x_size =
      naive_jacobian_matrix_.number_of_cameras() * params_per_camera +
      (naive_jacobian_matrix_.number_of_points() - number_of_gcps_) *
      params_per_point;
    
    if (i < y_feature_size)
    {
      return naive_jacobian_matrix_.coeff(i, j);
    }
    else if ((i - y_feature_size) == (j - naive_x_size))
    {
      return Scalar(1);
    }
    else
    {
      return Scalar(0);
    }
  }

private:
  NaiveJacobianMatrix naive_jacobian_matrix_;
  Index number_of_gcps_;
};

}
}
}

#endif