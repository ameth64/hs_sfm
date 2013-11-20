#ifndef _HS_SFM_SFM_FILE_IO_MATCHES_SAVER_HPP_
#define _HS_SFM_SFM_FILE_IO_MATCHES_SAVER_HPP_

#include <string>
#include <fstream>

#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{
namespace fileio
{

struct MatchesSaver
{
  typedef hs::sfm::MatchContainer Matches;
  typedef int Err;
  Err operator()(const std::string& matches_path,
                 const Matches& matches) const
  {
    std::ofstream matches_file(matches_path.c_str(), std::ios::out);
    if (!matches_file.is_open())
    {
      return -1;
    }

    Matches::const_iterator itr = matches.begin();
    Matches::const_iterator itr_end = matches.end();
    for (; itr != itr_end; ++itr)
    {
      matches_file<<(itr->first).first<<" "<<(itr->first).second<<"\n";
      size_t number_of_key_matches = itr->second.size();
      matches_file<<number_of_key_matches<<"\n";
      for (size_t i = 0; i < number_of_key_matches; i++)
      {
        matches_file<<itr->second[i].first<<" "<<
                      itr->second[i].second<<"\n";
      }
    }

    matches_file.close();

    return 0;
  }
};

}
}
}

#endif
