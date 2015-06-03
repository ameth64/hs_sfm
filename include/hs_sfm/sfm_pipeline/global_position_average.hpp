#ifndef _HS_SFM_SFM_PIPELINE_GLOBAL_POSITION_AVERAGE_HPP_
#define _HS_SFM_SFM_PIPELINE_GLOBAL_POSITION_AVERAGE_HPP_

#include <utility>

#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_pipeline/global_types.hpp"
#include "hs_sfm/sfm_pipeline/openmvg_utility.hpp"

namespace hs
{
namespace sfm
{
namespace pipeline
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
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::math::geometry::Rotation3D<Scalar> Rotation;
  typedef EIGEN_STD_MAP(size_t, Rotation) RotationContainer;

private:
  typedef hs::sfm::pipeline::OpenMVGUtility<Scalar> OpenMVGUtility;

public:
  int operator() (const KeysetContainer& keysets,
                  const hs::sfm::ObjectIndexMap& image_intrinsic_map,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const MatchContainer& matches,
                  const PositionPairContainer& position_pairs,
                  const RotationContainer& global_rotations,
                  IndexContainer& global_position_indices,
                  PositionContainer& global_positions) const
  {
    OpenMVGUtility openmvg_utility;

    //Fill sfm_data
    openMVG::SfM_Data sfm_data;
    openmvg_utility.FillSFMData(image_intrinsic_map, intrinsic_params_set,
                                sfm_data);

    //Fill normalized_features_provider
    openMVG::Features_Provider normalized_features_provider;
    openmvg_utility.FillNormalizedFeaturesProvider(
      keysets, image_intrinsic_map, intrinsic_params_set,
      normalized_features_provider);

    //Fill matches_provider
    openMVG::Matches_Provider matches_provider;
    openmvg_utility.FillMatchesProvider(matches, matches_provider);

    //Fill map_global_r
    openMVG::Hash_Map<openMVG::IndexT, openMVG::Mat3> map_global_r;
    openmvg_utility.FillMapGlobalR(global_rotations, map_global_r);

    openMVG::matching::PairWiseMatches triplet_wise_matches;
    openMVG::globalSfM::GlobalSfM_Translation_AveragingSolver solver;
    bool result = solver.Run(openMVG::globalSfM::TRANSLATION_AVERAGING_L1,
                             sfm_data,
                             &normalized_features_provider,
                             &matches_provider,
                             map_global_r,
                             triplet_wise_matches);

    global_position_indices.clear();
    global_positions.clear();
    auto itr_pose = sfm_data.poses.begin();
    auto itr_pose_end = sfm_data.poses.end();
    for (; itr_pose != itr_pose_end; ++itr_pose)
    {
      Position position;
      position << Scalar(itr_pose->second.center()[0]),
                  Scalar(itr_pose->second.center()[1]),
                  Scalar(itr_pose->second.center()[2]);
      global_position_indices.push_back(size_t(itr_pose->first));
      global_positions.push_back(position);
    }

    return result ? 0 : -1;
  }
};

}
}
}

#endif
