#ifndef _HS_SFM_SFM_PIPELINE_EPIPOLAR_EDGES_CALCULATOR_HPP_
#define _HS_SFM_SFM_PIPELINE_EPIPOLAR_EDGES_CALCULATOR_HPP_

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/essential/ematrix_5_points_calculator.hpp"
#include "hs_sfm/essential/ematrix_5_points_ransac_refiner.hpp"
#include "hs_sfm/essential/ematrix_extrinsic_params_points_calculator.hpp"
#include "hs_sfm/sfm_pipeline/global_types.hpp"

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
class EpipolarEdgesCalculator
{
public:
  typedef _Scalar Scalar;
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::sfm::pipeline::EpipolarEdge<Scalar> EpipolarEdge;
  typedef EIGEN_STD_VECTOR(EpipolarEdge) EpipolarEdgeContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

private:
  typedef typename EpipolarEdge::Index Index;
  typedef hs::sfm::essential::EMatrix5PointsRansacRefiner<Scalar>
          EMatrixRansacRefiner;
  typedef hs::sfm::essential::EMatrix5PointsCalculator<Scalar>
          EMatrixCalculator;
  typedef hs::sfm::essential::EMatrixExtrinsicParamsPointsCalculator<Scalar>
          ExtrinsicParamsPointsCalculator;
  typedef typename ExtrinsicParamsPointsCalculator::HKeyPair HKeyPair;
  typedef typename ExtrinsicParamsPointsCalculator::HKeyPairContainer
                   HKeyPairContainer;
  typedef typename ExtrinsicParamsPointsCalculator::EMatrix EMatrix;
  typedef typename IntrinsicParams::KMatrix KMatrix;
  typedef typename EMatrixRansacRefiner::IndexSet IndexSet;

public:
  int operator() (const ObjectIndexMap& image_intrinsic_map,
                  const hs::sfm::MatchContainer& matches,
                  const KeysetContainer& keysets,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  EpipolarEdgeContainer& epipolar_edges) const
  {
    EMatrixRansacRefiner ransac_refiner;
    ExtrinsicParamsPointsCalculator extrinsic_params_calculator;
    EMatrixCalculator calculator;
    auto itr_image_pair = matches.begin();
    auto itr_image_pair_end = matches.end();
    for (; itr_image_pair != itr_image_pair_end; ++itr_image_pair)
    {
      if (itr_image_pair->second.size() < 32) continue;
      HKeyPairContainer key_pairs;
      size_t image_id_left = itr_image_pair->first.first;
      size_t image_id_right = itr_image_pair->first.second;
      size_t intrinsic_id_left = image_intrinsic_map[image_id_left];
      size_t intrinsic_id_right = image_intrinsic_map[image_id_right];
      const IntrinsicParams& intrinsic_params_left =
        intrinsic_params_set[intrinsic_id_left];
      const IntrinsicParams& intrinsic_params_right =
        intrinsic_params_set[intrinsic_id_right];
      KMatrix K_left_inverse = intrinsic_params_left.GetKMatrix().inverse();
      KMatrix K_right_inverse = intrinsic_params_right.GetKMatrix().inverse();
      auto itr_key_pair = itr_image_pair->second.begin();
      auto itr_key_pair_end = itr_image_pair->second.end();
      for (; itr_key_pair != itr_key_pair_end; ++itr_key_pair)
      {
        size_t key_left_id = itr_key_pair->first;
        size_t key_right_id = itr_key_pair->second;
        HKeyPair key_pair;
        key_pair.first.segment(0, 2) = keysets[image_id_left][key_left_id];
        key_pair.first[2] = Scalar(1);
        key_pair.second.segment(0, 2) = keysets[image_id_right][key_right_id];
        key_pair.second[2] = Scalar(1);
        key_pair.first = K_left_inverse * key_pair.first;
        key_pair.second = K_right_inverse * key_pair.second;
        key_pairs.push_back(key_pair);
      }

      EpipolarEdge epipolar_edge;
      epipolar_edge.first_id = image_id_left;
      epipolar_edge.second_id = image_id_right;
      HKeyPairContainer key_pairs_refined;
      IndexSet inlier_indices;
      EMatrix e_matrix;
      if (ransac_refiner(key_pairs, 8/ intrinsic_params_left.focal_length(),
                         key_pairs_refined, inlier_indices,
                         e_matrix) != 0) continue;

      calculator(key_pairs_refined, e_matrix);

      PointContainer points_essential;
      if (extrinsic_params_calculator(
            e_matrix, key_pairs_refined,
            epipolar_edge.extrinsic_params_relative,
            points_essential) != 0) continue;

      epipolar_edges.push_back(epipolar_edge);
    }

    return 0;
  }
};

}
}
}

#endif
