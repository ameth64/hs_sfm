#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NOISED_Y_GENERATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NOISED_Y_GENERATOR_HPP_

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_feature_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BANaiveNoisedYGenerator
{
public:
  typedef int Err;
  typedef _Scalar Scalar;
  typedef BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::YVector YVector;
  typedef BANaiveFeatureCovarianceInverse<Scalar> YCovarianceInverse;

  typedef EIGEN_VECTOR(Scalar, 2) Vector2;

  Err operator()(const YVector& true_y,
                 const YCovarianceInverse& y_covariance_inverse,
                 YVector& noised_y) const
  {
    size_t number_of_features = true_y.rows() / 2;
    noised_y.resize(true_y.rows());
    for (size_t i = 0; i < number_of_features; i++)
    {
      Vector2 mean = true_y.segment(i * 2, 2);
      Vector2 noised_segment;
      hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
        mean, y_covariance_inverse.blocks[i].inverse(), noised_segment);
      noised_y.segment(i * 2, 2) = noised_segment;
    }

    return 0;
  }
};

}
}
}

#endif
