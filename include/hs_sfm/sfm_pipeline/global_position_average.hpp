#ifndef _HS_SFM_SFM_PIPELINE_GLOBAL_POSITION_AVERAGE_HPP_
#define _HS_SFM_SFM_PIPELINE_GLOBAL_POSITION_AVERAGE_HPP_

#include <utility>

#include "hs_sfm/sfm_pipeline/global_types.hpp"

namespace hs
{
namespace sfm
{
namespace global
{

template <typename _Scalar>
class GlobalPositionAverage
{
public:
  typedef _Scalar Scalar;
  typedef hs::sfm::pipeline::PositionPair<Scalar> PositionPair;
  typedef EIGEN_STD_VECTOR(PositionPair) PositionPairContainer;;
  typedef typename PositionPair::Position Position;
  typedef EIGEN_STD_VECTOR(Position) PositionContainer;
  typedef size_t Index;
  typedef std::vector<Index> IndexContainer;

public:
  int operator() (const PositionPairContainer& position_pairs,
                  IndexContainer& global_position_indices,
                  PositionContainer& global_positions) const
  {
    return 0;
  }
};

}
}
}

#endif
