#ifndef _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_NOISED_Y_GENERATOR_HPP_
#define _HS_SFM_PROJECTIVE_MLE_INTRINSIC_CONSTRAINED_NOISED_Y_GENERATOR_HPP_

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/projective/mle/intrinsic_constrained_vector_function.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_y_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

template <typename _Scalar>
class IntrinsicConstrainedNoisedYGenerator
{
public:
  typedef int Err;
  typedef _Scalar Scalar;
  typedef IntrinsicConstrainedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::YVector YVector;
  typedef IntrinsicConstrainedYCovarianceInverse<Scalar> YCovarianceInverse;

private:
  typedef typename VectorFunction::Index Index;
  typedef typename YCovarianceInverse::KeyBlock KeyBlock;

public:
  Err operator() (const YVector& true_y,
                  const YCovarianceInverse& y_covariance_inverse,
                  YVector& noised_y) const
  {
    Index number_of_keys = Index(y_covariance_inverse.key_blocks.size());
    noised_y.resize(number_of_keys * 2 + 5);
    for (Index i = 0; i < number_of_keys; i++)
    {
      const KeyBlock& key_block = y_covariance_inverse.key_blocks[i];
      EIGEN_VECTOR(Scalar, 2) key_noised;
      hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
        true_y.template segment<2>(i * 2), key_block,
        key_noised);
      noised_y.template segment<2>(i * 2) = key_noised;
    }

    hs::math::random::NormalRandomVar<Scalar, 1>::Generate(
      true_y[number_of_keys * 2 + 0],
      y_covariance_inverse.focal_length_stddev,
      noised_y[number_of_keys * 2 + 0]);

    hs::math::random::NormalRandomVar<Scalar, 1>::Generate(
      true_y[number_of_keys * 2 + 1],
      y_covariance_inverse.skew_stddev,
      noised_y[number_of_keys * 2 + 1]);

    hs::math::random::NormalRandomVar<Scalar, 1>::Generate(
      true_y[number_of_keys * 2 + 2],
      y_covariance_inverse.principal_point_x_stddev,
      noised_y[number_of_keys * 2 + 2]);

    hs::math::random::NormalRandomVar<Scalar, 1>::Generate(
      true_y[number_of_keys * 2 + 3],
      y_covariance_inverse.principal_point_y_stddev,
      noised_y[number_of_keys * 2 + 3]);

    hs::math::random::NormalRandomVar<Scalar, 1>::Generate(
      true_y[number_of_keys * 2 + 4],
      y_covariance_inverse.pixel_ratio_stddev,
      noised_y[number_of_keys * 2 + 4]);

    return 0;
  }
};

}
}
}

#endif