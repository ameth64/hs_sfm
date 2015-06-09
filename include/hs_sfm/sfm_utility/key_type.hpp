#ifndef _HS_SFM_UTILITY_KEY_TYPE_HPP_
#define _HS_SFM_UTILITY_KEY_TYPE_HPP_

#include <vector>

#include <cereal/types/vector.hpp>

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{

template <typename _Scalar>
class ImageKeys
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef _Scalar Scalar;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef size_t KeyId;
  typedef EIGEN_STD_VECTOR(Key) KeysContainer;

  ImageKeys() {}
  ImageKeys(const KeysContainer& keys) : keys_(keys) {}
  ImageKeys(size_t size) : keys_(size) {}

  void resize(size_t s) {keys_.resize(s);}
  size_t size() const {return keys_.size();}
  void clear() {keys_.swap(KeysContainer());}
  void AddKey(const Key& key) {keys_.push_back(key);}
  Key operator [] (KeyId i) const {return keys_[i];}
  Key& operator [] (KeyId i) {return keys_[i];}
  bool operator == (const ImageKeys& image_keys) const
  {
    return keys_ == image_keys.keys_;
  }

private:
  friend class cereal::access;

  template <class Archive>
  void serialize(Archive& archive)
  {
    archive(keys_);
  }

  KeysContainer keys_;
};

}
}

#endif
