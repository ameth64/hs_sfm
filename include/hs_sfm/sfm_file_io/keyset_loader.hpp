#ifndef _HS_SFM_SFM_FILE_IO_KEYSET_LOADER_HPP_
#define _HS_SFM_SFM_FILE_IO_KEYSET_LOADER_HPP_

#include <string>
#include <fstream>

#include "hs_sfm/sfm_utility/key_type.hpp"

namespace hs
{
namespace sfm
{
namespace fileio
{

template <typename _Scalar>
struct KeysetLoader
{
  typedef _Scalar Scalar;
  typedef ImageKeys<Scalar> Keys;
  typedef int Err;

  Err operator()(const std::string& keys_path,
                 Keys& keys) const
  {
    std::ifstream keys_file(keys_path.c_str(), std::ios::in);
    if (!keys_file.is_open())
    {
      return -1;
    }

    size_t number_of_keys;
    keys_file>>number_of_keys;
    keys.resize(number_of_keys);
    for (size_t i = 0; i < number_of_keys; i++)
    {
      keys_file>>keys[i][0]>>keys[i][1];
    }

    return 0;
  }
};

}
}
}

#endif
