#ifndef _HS_SFM_UTILITY_MATCHE_TYPE_HPP_
#define _HS_SFM_UTILITY_MATCHE_TYPE_HPP_

#include <vector>
#include <map>

namespace hs
{
namespace sfm
{

typedef std::vector<std::vector<int> > MatchGuide;
typedef std::pair<size_t, size_t> ImgPair;
typedef std::pair<size_t, size_t> KeyPair;
typedef std::vector<KeyPair> KeyPairContainer;
typedef std::map<ImgPair, KeyPairContainer> MatchContainer;

}
}

#endif
