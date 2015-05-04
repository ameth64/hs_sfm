#ifndef _HS_SFM_SFM_GLOBAL_GLOBAL_ROTATION_AVERAGE_HPP_
#define _HS_SFM_SFM_GLOBAL_GLOBAL_ROTATION_AVERAGE_HPP_

#include <utility>

#include "hs_sfm/sfm_pipeline/global_types.hpp"

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
class GlobalRotationAverage
{
public:
  typedef _Scalar Scalar;
  typedef hs::sfm::pipeline::RotationPair<Scalar> RotationPair;
  typedef EIGEN_STD_VECTOR(RotationPair) RotationPairContainer;
  typedef typename RotationPair::Rotation Rotation;
  typedef EIGEN_STD_VECTOR(Rotation) RotationContainer;
  typedef size_t Index;
  typedef std::vector<Index> IndexContainer;

public:
  int operator() (const RotationPairContainer& rotation_pairs,
                  IndexContainer& global_rotation_indices,
                  RotationContainer& global_rotation) const
  {
    return 0;
  }
};

}
}
}

#endif
