#ifndef _HS_SFM_SFM_FILE_IO_OBJECT_INDEX_LOADER_HPP_
#define _HS_SFM_SFM_FILE_IO_OBJECT_INDEX_LOADER_HPP_

#include <string>
#include <fstream>

#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{
namespace fileio
{

struct ObjectIndexLoader
{
  typedef int Err;
  typedef hs::sfm::ObjectIndexMap ObjectIndexMap;

  Err operator() (const std::string& object_index_path,
                  ObjectIndexMap& object_index_map) const
  {
    std::ifstream object_index_file(object_index_path.c_str(),
                                    std::ios::in);
    if (!object_index_file.is_open())
    {
      return -1;
    }

    size_t number_of_objects;
    object_index_file>>number_of_objects;
    object_index_map.Resize(number_of_objects);
    for (size_t i = 0; i < number_of_objects; i++)
    {
      size_t id;
      object_index_file>>id;
      if (id != 0)
      {
        object_index_map[i] = id - 1;
      }
    }

    return 0;
  }
};

}
}
}

#endif
