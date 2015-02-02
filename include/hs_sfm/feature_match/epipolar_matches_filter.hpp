#ifndef _HS_SFM_FEATURE_MATCH_EPIPOLAR_MATCHES_FILTER_HPP_
#define _HS_SFM_FEATURE_MATCH_EPIPOLAR_MATCHES_FILTER_HPP_

#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"

namespace hs
{
namespace sfm
{
namespace feature_match
{

template <typename _Scalar>
class EpipolarMatchesFilter
{
public:
  typedef _Scalar Scalar;
  typedef ImageKeys<Scalar> ImageKeysType;
  typedef EIGEN_STD_VECTOR(ImageKeysType) ImageKeysContainer;

public:
  EpipolarMatchesFilter(Scalar filter_threshold);
    : filter_threshold_(filter_threshold) {}

  int operator() (const ImageKeysContainer& image_keys_set,
                  const MatchContainer& matches,
                  MatchContainer& matches_filtered) const
  {
    return 0;
  }

private:
  Scalar filter_threshold_;
};

}
}
}

#endif
