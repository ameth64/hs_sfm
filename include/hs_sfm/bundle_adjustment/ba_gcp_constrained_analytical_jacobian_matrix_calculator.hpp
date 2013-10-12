#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_ANALYTICAL_JACOBIAN_MATRIX_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_ANALYTICAL_JACOBIAN_MATRIX_CALCULATOR_HPP_

#include "hs_sfm/bundle_adjustment/ba_naive_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_jacobian_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VectorFunction>
class BAGCPConstrainedAnalyticalJacobianMatrixCalculator;

template <typename _Scalar>
class BAGCPConstrainedAnalyticalJacobianMatrixCalculator<
  BAGCPConstrainedVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef BAGCPConstrainedVectorFunction<Scalar> VectorFunction;
  typedef BANaiveVectorFunction<Scalar> NaiveVectorFunction;
  typedef BANaiveAnalyticalJacobianMatrixCalculator<NaiveVectorFunction>
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