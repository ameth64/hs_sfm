#ifndef _HS_SFM_FEATURE_MATCH_FLANN_MATCHER_HPP_
#define _HS_SFM_FEATURE_MATCH_FLANN_MATCHER_HPP_

namespace hs
{
namespace sfm
{
namespace feature_match
{

class HS_EXPORT FLANNMatcher
{
public:
  FLANNMatcher(float match_threshold = 0.6f);
};

}
}
}

#endif
