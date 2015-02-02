#ifndef _HS_SFM_INCREMENTAL_BEST_PAIR_SELECTOR_HPP_
#define _HS_SFM_INCREMENTAL_BEST_PAIR_SELECTOR_HPP_

#include <map>
#include <limits>

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/triangulate/multiple_view_vector_function.hpp"
#include "hs_sfm/triangulate/multiple_view_levenberg_marquardt_optimizor.hpp"
#include "hs_sfm/essential/ematrix_5_points_calculator.hpp"
#include "hs_sfm/essential/ematrix_5_points_ransac_refiner.hpp"
#include "hs_sfm/essential/ematrix_extrinsic_params_points_calculator.hpp"
#if 1
#include "hs_sfm/sfm_utility/projective_functions.hpp"
#endif

namespace hs
{
namespace sfm
{
namespace incremental
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
  typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;
  typedef typename EMatrixCalculator::EMatrix EMatrix;
  typedef hs::sfm::essential::EMatrixExtrinsicParamsPointsCalculator<Scalar>
          ExtrinsicParamsPointsCalculator;
#if 1
  typedef ProjectiveFunctions<Scalar> ProjectiveFunctionsType;
  typedef typename ProjectiveFunctionsType::Key Key;
#endif

public:
  BestPairSelector(size_t min_number_of_pair_matches)
    : min_number_of_pair_matches_(min_number_of_pair_matches) {}

  Err operator()(const ImageKeysetContainer& image_keysets,
                 const MatchContainer& matches,
                 const IntrinsicParamsContainer& intrinsic_params_set,
                 const ObjectIndexMap& image_intrinsic_map,
                 size_t& best_identity_id,
                 size_t& best_relative_id,
                 ExtrinsicParams& relative_extrinsic_params,
                 PointContainer& points) const
  {
    auto image_pair_itr = matches.begin();
    auto image_pair_itr_end = matches.end();
    auto image_pair_itr_best = matches.begin();
    size_t max_number_of_pair_matches = 0;
    for (; image_pair_itr != image_pair_itr_end; ++image_pair_itr)
    {
      if (image_pair_itr->second.size() > max_number_of_pair_matches)
      {
        image_pair_itr_best = image_pair_itr;
        max_number_of_pair_matches = image_pair_itr->second.size();
      }
    }
    if (max_number_of_pair_matches < min_number_of_pair_matches_)
    {
      return -1;
    }
    if (image_pair_itr_best == image_pair_itr_end)
    {
      return -1;
    }
    image_pair_itr_best = matches.find(ImagePair(12, 13));
    std::cout<<"image_pair_itr_best:"<<image_pair_itr_best->first.first<<" "
                                     <<image_pair_itr_best->first.second<<"\n";
    auto key_pair_itr = image_pair_itr_best->second.begin();
    auto key_pair_itr_end = image_pair_itr_best->second.end();
    HKeyPairContainer key_pairs;
    size_t image_left_id = image_pair_itr_best->first.first;
    size_t intrinsic_id_left = image_intrinsic_map[image_left_id];
    size_t image_right_id = image_pair_itr_best->first.second;
    size_t intrinsic_id_right = image_intrinsic_map[image_right_id];
    const IntrinsicParams& intrinsic_params_left =
      intrinsic_params_set[intrinsic_id_left];
    const IntrinsicParams& intrinsic_params_right =
      intrinsic_params_set[intrinsic_id_right];
    KMatrix K_left_inverse = intrinsic_params_left.GetKMatrix().inverse();
    KMatrix K_right_inverse = intrinsic_params_right.GetKMatrix().inverse();
    for (; key_pair_itr != key_pair_itr_end; ++key_pair_itr)
    {
      size_t key_left_id = key_pair_itr->first;
      size_t key_right_id = key_pair_itr->second;
      HKeyPair key_pair;
      key_pair.first.segment(0, 2) =
        image_keysets[image_left_id][key_left_id];
      key_pair.first[2] = Scalar(1);
      key_pair.second.segment(0, 2) =
        image_keysets[image_right_id][key_right_id];
      key_pair.second[2] = Scalar(1);
      key_pair.first = K_left_inverse * key_pair.first;
      key_pair.second = K_right_inverse * key_pair.second;
      key_pairs.push_back(key_pair);
    }
    EMatrix e_matrix;
    EMatrixCalculator ematrix_calculator;
    if (ematrix_calculator(key_pairs, e_matrix) != 0)
    {
      return -1;
    }

    //通过E矩阵计算影像对的相对外方位元素
    ExtrinsicParamsPointsCalculator extrinsic_points_calculator;
    if (extrinsic_points_calculator(e_matrix,
                                    key_pairs,
                                    relative_extrinsic_params,
                                    points) != 0)
    {
      return -1;
    }
    best_identity_id = image_left_id;
    best_relative_id = image_right_id;

//    auto image_pair_itr = matches.begin();
//    auto image_pair_itr_end = matches.end();
//    EMatrixCalculator ematrix_calculator;
//    Scalar min_mean_height = -std::numeric_limits<Scalar>::max();
//    for (; image_pair_itr != image_pair_itr_end; image_pair_itr++)
//    {
//      if (image_pair_itr->second.size() > min_number_of_pair_matches_)
//      {
//        auto key_pair_itr = image_pair_itr->second.begin();
//        auto key_pair_itr_end = image_pair_itr->second.end();
//        HKeyPairContainer key_pairs;
//        size_t image_left_id = image_pair_itr->first.first;
//        size_t intrinsic_id_left = image_intrinsic_map[image_left_id];
//        size_t image_right_id = image_pair_itr->first.second;
//        size_t intrinsic_id_right = image_intrinsic_map[image_right_id];
//        const IntrinsicParams& intrinsic_params_left =
//          intrinsic_params_set[intrinsic_id_left];
//        const IntrinsicParams& intrinsic_params_right =
//          intrinsic_params_set[intrinsic_id_right];
//        KMatrix K_left_inverse = intrinsic_params_left.GetKMatrix().inverse();
//        KMatrix K_right_inverse = intrinsic_params_right.GetKMatrix().inverse();
//        for (; key_pair_itr != key_pair_itr_end; ++key_pair_itr)
//        {
//          size_t key_left_id = key_pair_itr->first;
//          size_t key_right_id = key_pair_itr->second;
//          HKeyPair key_pair;
//          key_pair.first.segment(0, 2) =
//            image_keysets[image_left_id][key_left_id];
//          key_pair.first[2] = Scalar(1);
//          key_pair.second.segment(0, 2) =
//            image_keysets[image_right_id][key_right_id];
//          key_pair.second[2] = Scalar(1);
//          key_pair.first = K_left_inverse * key_pair.first;
//          key_pair.second = K_right_inverse * key_pair.second;
//          key_pairs.push_back(key_pair);
//        }
//
//        EMatrix e_matrix;
//        if (ematrix_calculator(key_pairs, e_matrix) != 0)
//        {
//          continue;
//        }
//
//        //通过E矩阵计算影像对的相对外方位元素
//        ExtrinsicParamsPointsCalculator extrinsic_points_calculator;
//        ExtrinsicParams extrinsic_params_pair;
//        PointContainer points_pair;
//        if (extrinsic_points_calculator(e_matrix,
//                                        key_pairs,
//                                        extrinsic_params_pair,
//                                        points_pair) != 0)
//        {
//          continue;
//        }
//
//        //计算平均基高比
//        Scalar mean_height = Scalar(0);
//        size_t number_of_points = points_pair.size();
//        for (size_t i = 0; i < number_of_points; i++)
//        {
//          mean_height += std::abs(points_pair[i][2]);
//#if 1
//          ExtrinsicParams extrinsic_params_left;
//          Key key_left_projected =
//            ProjectiveFunctionsType::WorldPointProjectToImageKey(
//              intrinsic_params_left, extrinsic_params_left, points_pair[i]);
//          Key key_right_projected =
//            ProjectiveFunctionsType::WorldPointProjectToImageKey(
//              intrinsic_params_right, extrinsic_params_pair, points_pair[i]);
//          size_t key_left_id = image_pair_itr->second[i].first;
//          size_t key_right_id = image_pair_itr->second[i].second;
//          Key key_left_predicated =
//            image_keysets[image_left_id][key_left_id];
//          Key key_right_predicated =
//            image_keysets[image_right_id][key_right_id];
//
//          Key diff_left = key_left_projected - key_left_predicated;
//          Key diff_right = key_right_projected - key_right_predicated;
//          //std::cout<<"diff_left:\n"<<diff_left<<"\n";
//          //std::cout<<"diff_right:\n"<<diff_right<<"\n";
//#endif
//        }
//        mean_height /= Scalar(number_of_points);
//
//        if (mean_height > min_mean_height)
//        {
//          best_identity_id = image_left_id;
//          best_relative_id = image_right_id;
//          relative_extrinsic_params = extrinsic_params_pair;
//          points.swap(points_pair);
//          min_mean_height = mean_height;
//        }
//      }
//    }

    return 0;
  }

private:
  size_t min_number_of_pair_matches_;
};

}
}
}

#endif
