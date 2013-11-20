#ifndef _HS_SFM_SFM_FILE_IO_KEYSET_SAVER_HPP_
#define _HS_SFM_SFM_FILE_IO_KEYSET_SAVER_HPP_

#include <string>
#include <fstream>
#include <iomanip>

#include "hs_sfm/sfm_utility/key_type.hpp"

namespace hs
{
namespace sfm
{
namespace fileio
{

template <typename _Scalar>
struct KeysetSaver
{
  typedef _Scalar Scalar;
  typedef ImageKeys<Scalar> Keys;
  typedef int Err;

  Err operator()(const std::string& keys_path, 
                 const Keys& keys) const
  {
    std::ofstream keys_file(keys_path.c_str(), std::ios::out);
    if (!keys_file.is_open())
    {
      return -1;
    }

    keys_file.setf(std::ios::fixed);
    keys_file<<std::setprecision(6);

    size_t number_of_keys = keys.size();
    keys_file<<number_of_keys<<'\n';
    for (size_t i = 0; i < number_of_keys; i++)
    {
      keys_file<<keys[i][0]<<' '<<keys[i][1]<<'\n';
    }

    return 0;
  }
};

}
}
}

#endif
