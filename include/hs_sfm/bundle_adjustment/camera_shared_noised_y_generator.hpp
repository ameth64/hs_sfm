#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NOISED_Y_GENERATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NOISED_Y_GENERATOR_HPP_

#include <cmath>

#include "hs_math/random/normal_random_var.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_y_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedNoisedYGenerator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::YVector YVector;
  typedef CameraSharedYCovarianceInverse<Scalar> YCovarianceInverse;

private:
  typedef typename VectorFunction::Index Index;
  typedef EIGEN_VECTOR(Scalar, VectorFunction::params_per_key_) Key;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) Constraint;
  typedef typename YCovarianceInverse::KeyBlock KeyBlock;
  typedef hs::math::random::NormalRandomVar<Scalar,
                                            VectorFunction::params_per_key_>
          KeyGenerator;
  typedef hs::math::random::NormalRandomVar<Scalar, 1>
          ConstraintGenerator;
public:
  Err operator() (const YVector& true_y,
                  const YCovarianceInverse& y_covariance_inverse,
                  YVector& noised_y) const
  {
    Index y_size = y_covariance_inverse.GetYSize();
    if (y_size != true_y.rows())
    {
      return -1;
    }
    noised_y.resize(y_size);
    size_t number_of_keys = y_covariance_inverse.NumberOfKeys();

    for (size_t i = 0; i < number_of_keys; i++)
    {
      Key true_key = true_y.segment(i * VectorFunction::params_per_key_,
                                    VectorFunction::params_per_key_);
      KeyBlock key_covariance = y_covariance_inverse.GetKeyBlock(i).inverse();
      Key noised_key;
      KeyGenerator::Generate(true_key, key_covariance, noised_key);
      noised_y.segment(i * VectorFunction::params_per_key_,
                       VectorFunction::params_per_key_) = noised_key;
    }

    Index constraints_size = y_covariance_inverse.GetConstraintsSize();
    Index key_params_size = y_covariance_inverse.GetKeyParamsSize();
    for (Index constraint_id = 0; constraint_id < constraints_size;
         constraint_id++)
    {
      Scalar sigma =
        1 / std::sqrt(y_covariance_inverse.GetConstraint(constraint_id));
      Scalar mean = true_y[key_params_size + constraint_id];
      ConstraintGenerator::Generate(
        mean, sigma, noised_y[key_params_size + constraint_id]);
    }

    return 0;
  }
};

}
}
}

#endif
