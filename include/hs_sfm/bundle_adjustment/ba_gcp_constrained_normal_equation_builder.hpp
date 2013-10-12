#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_NORMAL_MATRIX_BUILDER_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_NORMAL_MATRIX_BUILDER_HPP_

#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_jacobian_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_y_covariance_inverse.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_equation_builder.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_gradient.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BAGCPConstrainedNormalEquationBuilder
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef BAGCPConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;
  typedef BAGCPConstrainedJacobianMatrix<Scalar, Index,
                                         VectorFunction::params_per_feature_,
                                         VectorFunction::params_per_camera_,
                                         VectorFunction::params_per_point_>
          JacobianMatrix;
  typedef BAGCPConstrainedYCovarianceInverse<Scalar> YCovarianceInverse;
  typedef BANaiveNormalEquationBuilder<Scalar> NaiveBuilder;
  typedef typename NaiveBuilder::Residuals Residuals;
  typedef typename NaiveBuilder::NormalMatrix NormalMatrix;
  typedef typename NaiveBuilder::Gradient Gradient;

  typedef typename NormalMatrix::PointBlock PointBlock;
  typedef typename Gradient::PointSegment PointSegment;

  Err operator() (const JacobianMatrix& jacobian_matrix,
                  const Residuals& residuals,
                  const YCovarianceInverse& y_covariance_inverse,
                  NormalMatrix& normal_matrix,
                  Gradient& gradient) const
  {
    NaiveBuilder naive_builder;
    Err result = naive_builder(jacobian_matrix.naive_jacobian_matrix(),
                               residuals,
                               y_covariance_inverse.naive_y_covariance_inverse,
                               normal_matrix,
                               gradient);
    if (result != 0)
    {
      return result;
    }

    Index gcp_id_start = normal_matrix.number_of_points - 
                         jacobian_matrix.number_of_gcps();
    Index naive_y_size = 
      y_covariance_inverse.naive_y_covariance_inverse.blocks.size() *
      VectorFunction::params_per_feature_;
    for (Index i = 0; i < jacobian_matrix.number_of_gcps(); i++)
    {
      normal_matrix.point_blocks[gcp_id_start + i] +=
        y_covariance_inverse.gcp_blocks[i];
      gradient.point_segments[gcp_id_start + i] +=
        y_covariance_inverse.gcp_blocks[i] *
        residuals.segment(naive_y_size + i * VectorFunction::params_per_point_,
                          VectorFunction::params_per_point_);
    }

    return 0;
  }
};

}
}
}

#endif