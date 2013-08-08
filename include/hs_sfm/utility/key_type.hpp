#ifndef _HS_SFM_UTILITY_KEY_TYPE_HPP_
#define _HS_SFM_UTILITY_KEY_TYPE_HPP_

#include <vector>

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
  typedef EIGEN_VEC(Scalar, 2) Key;
  typedef size_t KeyId;
  typedef EIGEN_VECTOR(Key) KeysContainer;

  ImageKeys() {}
  ImageKeys(const KeysContainer& keys) : m_keys(keys) {}
  ImageKeys(size_t size) : m_keys(size) {}

  void resize(size_t s) {m_keys.resize(s);}
  size_t size() const {return m_keys.size();}
  void clear() {m_keys.swap(KeysContainer());}
  Key operator [] (KeyId i) const {return m_keys[i];}
  Key& operator [] (KeyId i) {return m_keys[i];}
  bool operator == (const ImageKeys& imgKeys)
  {
    return m_keys == imgKeys.m_keys;
  }

private:
  KeysContainer m_keys;
};

}
}

#endif
