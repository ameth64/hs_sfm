#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_FORWARD_FINITE_DIFFERENCE_JACOBIAN_MATRIX_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_FORWARD_FINITE_DIFFERENCE_JACOBIAN_MATRIX_CALCULATOR_HPP_

#include "hs_sfm/bundle_adjustment/ba_naive_forward_finite_difference_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_jacobian_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{
template <typename _VectorFunction>
class BAGCPConstrainedForwardFiniteDifferenceJacobianMatrixCalculator;

/**
 *  Forward Finite Difference计算带GCP约束的Bundle Adjustment函数的Jacobian矩阵
 */
template <typename _Scalar>
class BAGCPConstrainedForwardFiniteDifferenceJacobianMatrixCalculator<
  BAGCPConstrainedVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef BAGCPConstrainedVectorFunction<Scalar> VectorFunction;
  typedef BANaiveVectorFunction<Scalar> NaiveVectorFunction;
  typedef BANaiveForwardFiniteDifferenceJacobianMatrixCalculator<
            NaiveVectorFunction>
          NaiveJacobianMatrixCalculator;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;
  typedef BAGCPConstrainedJacobianMatrix<Scalar, Index,
                                         VectorFunction::params_per_feature_,
                                         VectorFunction::params_per_camera_,
                                         VectorFunction::params_per_point_>
          JacobianMatrix;
public:
  BAGCPConstrainedForwardFiniteDifferenceJacobianMatrixCalculator() {}
  BAGCPConstrainedForwardFiniteDifferenceJacobianMatrixCalculator(
    Scalar delta, Scalar min_delta)
    : naive_jacobian_matrix_calculator_(delta, min_delta) {}

  Err operator() (const VectorFunction& vector_function,
                  const XVector& x,
                  JacobianMatrix& jacobian_matrix) const
  {
    jacobian_matrix.set_number_of_gcps(vector_function.number_of_gcps());
    return naive_jacobian_matrix_calculator_(
      NaiveVectorFunction(vector_function),
      x, jacobian_matrix.naive_jacobian_matrix());
  }
private:
  NaiveJacobianMatrixCalculator naive_jacobian_matrix_calculator_;
};

}
}
}

#endif