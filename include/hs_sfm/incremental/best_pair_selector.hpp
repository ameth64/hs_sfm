#ifndef _HS_SFM_INCREMENTAL_BEST_PAIR_SELECTOR_HPP_
#define _HS_SFM_INCREMENTAL_BEST_PAIR_SELECTOR_HPP_

#include <map>

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/triangulate/multiple_view_vector_function.hpp"
#include "hs_sfm/triangulate/multiple_view_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/essential/ematrix_5_points_calculator.hpp"
#include "hs_sfm/essential/ematrix_5_points_ransac_refiner.hpp"

namespace hs
{
namespace sfm
{
namespace projective
{

/**
 *  选择最优的两张照片。
 *
 *  选择满足最小匹配数，并且基高比最大的两张照片。
 */
template <typename _Scalar>
class BestPairSelector
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef ImageKeys<Scalar> ImageKeyset;
  typedef EIGEN_STD_VECTOR(ImageKeyset) ImageKeysetContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

private:
  typedef std::map<ImagePair, Scalar> ScoreContainer;
  typedef hs::sfm::essential::EMatrix5PointsCalculator<Scalar>
          EMatrixCalculator;
  typedef typename EMatrixCalculator::HKey HKey;
  typedef typename EMatrixCalculator::HKeyPair HKeyPair;
  typedef typename EMatrixCalculator::HKeyPairContainer HKeyPairContainer;
  typedef typename IntrinsicParams::KMatrix KMatrix;
  typedef typename EIGEN_MATRIX(Scalar, 3, 3) RMatrix;

public:
  BestPairSelector(size_t min_number_of_pair_matches)
    : min_number_of_pair_matches_(min_number_of_pair_matches) {}

  Err operator()(const ImageKeysetContainer& image_keysets,
                 const MatchContainer& matches,
                 const IntrinsicParamsContainer& intrinsic_params_set,
                 size_t& best_identity_id,
                 size_t& best_relative_id,
                 ExtrinsicParams& relative_extrinsic_params,
                 PointContainer& points) const
  {
    auto image_pair_itr = matches.begin();
    auto image_pair_itr_end = matches.end();
    EMatrixCalculator ematrix_calculator;
    for (; image_pair_itr != image_pair_itr_end; image_pair_itr++)
    {
      if (image_pair_itr->second.size() > min_number_of_pair_matches_)
      {
        auto key_pair_itr = image_pair_itr->second.begin();
        auto key_pair_itr_end = image_pair_itr->second.end();
        HKeyPairContainer key_pairs;
        size_t image_left_id = key_pair_itr->first.first;
        size_t image_right_id = key_pair_itr->first.second;
        const IntrinsicParams& intrinsic_params_left =
          intrinsic_params_set[image_left_id];
        KMatrix K_left_inverse = intrinsic_params_left.GetKMatrix().inverse();
        KMatrix K_right_inverse = intrinsic_params_right.GetKMatrix().inverse();
        const IntrinsicParams& intrinsic_params_right =
          intrinsic_params_set[image_right_id];
        for (; key_pair_itr != key_pair_itr_end; ++key_pair_itr)
        {
          size_t key_left_id = key_pair_itr->first;
          size_t key_right_id = key_pair_itr->second;
          HKeyPair key_pair;
          key_pair.first = image_keysets[image_left_id][key_left_id];
          key_pair.second = image_keysets[image_right_id][key_right_id];
          key_pair.first = K_left_inverse * key_pair.first;
          key_pair.second = K_right_inverse * key_pair.second;
          key_pairs.push_back(key_pair);
        }

        EMatrix e_matrix;
        if (ematrix_calculator(key_pairs, e_matrix) != 0)
        {
          return -1;
        }

        //通过E矩阵计算影像对的相对外方位元素
      }
    }

    return 0;
  }

private:
  size_t min_number_of_pair_matches_;
};

}
}
}

#endif