#ifndef _HS_SFM_SFM_FILE_IO_MATCHES_LOADER_HPP_
#define _HS_SFM_SFM_FILE_IO_MATCHES_LOADER_HPP_

#include <string>
#include <fstream>

#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{
namespace fileio
{

struct MatchesLoader
{
  typedef hs::sfm::ImagePair ImagePair;
  typedef hs::sfm::KeyPair KeyPair;
  typedef hs::sfm::KeyPairContainer KeyPairContainer;
  typedef hs::sfm::MatchContainer MatchContainer;
  typedef int Err;
  Err operator()(const std::string& matches_path,
                 MatchContainer& matches) const
  {
    std::ifstream matches_file(matches_path.c_str(), std::ios::in);
    if (!matches_file.is_open())
    {
      return -1;
    }

    matches.clear();
    while (!matches_file.eof())
    {
      char buf[512];
      matches_file.getline(buf, 512);
      if (buf[0] == 0)
      {
        break;
      }
      std::stringstream ss(buf);
      ImagePair image_pair;
      ss>>image_pair.first>>image_pair.second;
      size_t number_of_image_pairs;
      matches_file.getline(buf, 512);
      ss.clear();
      ss.str(buf);
      ss>>number_of_image_pairs;
      KeyPairContainer key_pairs(number_of_image_pairs);
      for (size_t i = 0; i < number_of_image_pairs; i++)
      {
        matches_file.getline(buf, 512);
        ss.clear();
        ss.str(buf);
        ss>>key_pairs[i].first>>key_pairs[i].second;
      }
      matches[image_pair] = key_pairs;
    }

    return 0;
  }
};

}
}
}

#endif
