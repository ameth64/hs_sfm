#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_NOISED_Y_GENERATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_NOISED_Y_GENERATOR_HPP_

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_gcp_constrained_y_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BAGCPConstrainedNoisedYGenerator
{
public:
  typedef int Err;
  typedef _Scalar Scalar;
  typedef BAGCPConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::YVector YVector;
  typedef BAGCPConstrainedYCovarianceInverse<Scalar> YCovarianceInverse;

private:
  typedef EIGEN_VECTOR(Scalar, VectorFunction::params_per_feature_) Feature;
  typedef EIGEN_VECTOR(Scalar, VectorFunction::params_per_point_) Point;

public:
  Err operator() (const YVector& true_y,
                  const YCovarianceInverse& y_covariance_inverse,
                  YVector& noised_y) const
  {
    size_t number_of_features =
      y_covariance_inverse.naive_y_covariance_inverse.blocks.size();
    size_t number_of_gcps =
      y_covariance_inverse.gcp_blocks.size();

    Index naive_y_size =
      number_of_features * VectorFunction::params_per_feature_;
    Index y_size = naive_y_size +
                   number_of_gcps * VectorFunction::params_per_point_;
    noised_y.resize(y_size);

    for (size_t i = 0; i < number_of_features; i++)
    {
      Feature mean = true_y.segment(i * VectorFunction::params_per_feature_,
                                    VectorFunction::params_per_feature_);
      Feature noised_feature;
      hs::math::random::NormalRandomVar<
        Scalar, VectorFunction::params_per_feature_>::Generate(
        mean,
        y_covariance_inverse.naive_y_covariance_inverse.blocks[i].inverse(),
        noised_feature);
      noised_y.segment(i * VectorFunction::params_per_feature_,
                       VectorFunction::params_per_feature_) = noised_feature;
    }

    for (size_t i = 0; i < number_of_gcps; i++)
    {
      Index gcp_id_start = naive_y_size + 
                           i * VectorFunction::params_per_point_;
      Point mean = true_y.segment(gcp_id_start,
                                  VectorFunction::params_per_point_);
      Point noised_point;
      hs::math::random::NormalRandomVar<
        Scalar, VectorFunction::params_per_point_>::Generate(
        mean, y_covariance_inverse.gcp_blocks[i].inverse(),
        noised_point);
      noised_y.segment(gcp_id_start, VectorFunction::params_per_point_) =
        noised_point;
    }

    return 0;
  }
};

}
}
}

#endif