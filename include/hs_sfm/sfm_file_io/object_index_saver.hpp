#ifndef _HS_SFM_SFM_FILE_IO_OBJECT_INDEX_SAVER_HPP_
#define _HS_SFM_SFM_FILE_IO_OBJECT_INDEX_SAVER_HPP_

#include <string>
#include <fstream>

#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{
namespace fileio
{

struct ObjectIndexSaver
{
  typedef int Err;
  typedef hs::sfm::ObjectIndexMap ObjectIndexMap;

  Err operator() (const std::string& object_index_path,
                  const ObjectIndexMap& object_index_map) const
  {
    std::ofstream object_index_file(object_index_path.c_str(),
                                    std::ios::out);
    if (!object_index_file.is_open())
    {
      return -1;
    }

    size_t number_of_objects = object_index_map.Size();
    object_index_file<<number_of_objects<<"\n";
    for (size_t i = 0; i < number_of_objects; i++)
    {
      if (object_index_map.IsValid(i))
      {
        object_index_file<<object_index_map[i] + 1<<"\n";
      }
      else
      {
        object_index_file<<"0\n";
      }
    }

    return 0;
  }
};

}
}
}

#endif
