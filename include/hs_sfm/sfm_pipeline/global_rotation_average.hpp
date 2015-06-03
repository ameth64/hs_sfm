#ifndef _HS_SFM_SFM_GLOBAL_GLOBAL_ROTATION_AVERAGE_HPP_
#define _HS_SFM_SFM_GLOBAL_GLOBAL_ROTATION_AVERAGE_HPP_

#include <utility>

#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"

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
                  RotationContainer& global_rotations) const
  {
    typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;
    openMVG::rotation_averaging::RelativeRotations relative_rotations(
      rotation_pairs.size());

    for (size_t i = 0; i < rotation_pairs.size(); i++)
    {
      relative_rotations[i].i = openMVG::IndexT(rotation_pairs[i].first_id);
      relative_rotations[i].j = openMVG::IndexT(rotation_pairs[i].second_id);
      relative_rotations[i].weight = float(rotation_pairs[i].weight);
      RMatrix r = rotation_pairs[i].rotation;
      relative_rotations[i].Rij = openMVG::Mat3(r);
    }

    openMVG::globalSfM::ERelativeRotationInferenceMethod inference_method =
      openMVG::globalSfM::TRIPLET_ROTATION_INFERENCE_COMPOSITION_ERROR;
    openMVG::globalSfM::GlobalSfM_Rotation_AveragingSolver solver;
    openMVG::Hash_Map<openMVG::IndexT, openMVG::Mat3> map_global_r;
    const bool result = solver.Run(openMVG::globalSfM::ROTATION_AVERAGING_L1,
                                   inference_method,
                                   relative_rotations,
                                   map_global_r);

    global_rotation_indices.clear();
    global_rotations.clear();
    auto itr_global_r = map_global_r.begin();
    auto itr_global_r_end = map_global_r.end();
    for (; itr_global_r != itr_global_r_end; ++itr_global_r)
    {
      global_rotation_indices.push_back(size_t(itr_global_r->first));
      RMatrix r = itr_global_r->second;
      global_rotations.push_back(Rotation(r));
    }

    return result ? 0 : -1;
  }
};

}
}
}

#endif
