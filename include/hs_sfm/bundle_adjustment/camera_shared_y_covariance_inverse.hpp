#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_Y_COVARIANCE_INVERSE_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_Y_COVARIANCE_INVERSE_HPP_

#include <utility>
#include <limits>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_common_types.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedYCovarianceInverse
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef int Err;
  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;

  typedef EIGEN_MATRIX(Scalar, VectorFunction::params_per_key_,
                               VectorFunction::params_per_key_)
          KeyBlock;
  typedef EIGEN_STD_VECTOR(KeyBlock) KeyBlockContainer;

private:
  typedef std::vector<Scalar> ConstraintSegment;

public:
  void Clear()
  {
    KeyBlockContainer().swap(key_blocks_);
    constraints_.resize(0);
  }

  Err SetKeysUniformStdDev(Scalar key_std_dev,
                           size_t number_of_keys)
  {
    if (key_std_dev == Scalar(0)) return -1;
    return 0;
  }

  Err SetKeysUniformCovariance(const KeyBlock& key_covariance,
                               size_t number_of_keys)
  {
    KeyBlock key_covariance_inverse = key_covariance.inverse();
    key_blocks_.resize(number_of_keys, key_covariance_inverse);
    return 0;
  }

  Err AddKeyCovariance(const KeyBlock& key_covariance)
  {
    KeyBlock key_covariance_inverse = key_covariance.inverse();
    key_blocks_.push_back(key_covariance_inverse);
    return 0;
  }

  const KeyBlock& GetKeyBlock(size_t key_id) const
  {
    return key_blocks_[key_id];
  }
  KeyBlock& GetKeyBlock(size_t key_id)
  {
    return key_blocks_[key_id];
  }

  void AddConstraint(Scalar constraint)
  {
    constraints_.push_back(constraint);
  }

  Scalar GetConstraint(Index constaint_id) const
  {
    return constraints_[constaint_id];
  }

  Scalar coeff(Index i, Index j) const
  {
    Index key_params_size =
      Index(key_blocks_.size()) * VectorFunction::params_per_key_;
    Index y_size = GetYSize();
    if (i < 0 || i >= y_size || j < 0 || j >= y_size)
    {
      return std::numeric_limits<Scalar>::signaling_NaN();
    }
    if (i < key_params_size)
    {
      //keys block
      Index block_id = i / VectorFunction::params_per_key_;
      Index row_id = i % VectorFunction::params_per_key_;
      Index block_begin = block_id * VectorFunction::params_per_key_;
      Index block_end = block_begin + VectorFunction::params_per_key_;
      if (j >= block_begin && j < block_end)
      {
        Index col_id = j - block_begin;
        return key_blocks_[block_id](row_id, col_id);
      }
      else
      {
        return Scalar(0);
      }
    }
    else
    {
      if (i == j)
      {
        Index offset = (i - key_params_size);
        return constraints_[offset];
      }
      else
      {
        return Scalar(0);
      }
    }
  }

  Index GetKeyParamsSize() const
  {
    return Index(key_blocks_.size()) * VectorFunction::params_per_key_;
  }

  Index GetConstraintsSize() const
  {
    return Index(constraints_.size());
  }

  Index GetYSize() const
  {
    return GetKeyParamsSize() + GetConstraintsSize();
  }

  size_t NumberOfKeys() const
  {
    return key_blocks_.size();
  }

private:
  KeyBlockContainer key_blocks_;
  ConstraintSegment constraints_;
};

}
}
}

#endif
