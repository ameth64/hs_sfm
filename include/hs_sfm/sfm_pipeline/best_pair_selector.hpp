#ifndef _HS_SFM_SFM_PIPELINE_BEST_PAIR_SELECTOR_HPP_
#define _HS_SFM_SFM_PIPELINE_BEST_PAIR_SELECTOR_HPP_

#include <map>
#include <limits>
#include <algorithm>

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
namespace pipeline
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
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;	//相机外参数矩阵
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;	//相机内参数矩阵
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;	//内参数矩阵容器
  typedef ImageKeys<Scalar> ImageKeyset;
  typedef EIGEN_STD_VECTOR(ImageKeyset) ImageKeysetContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

private:
  typedef std::map<ImagePair, Scalar> ScoreContainer;	//对每个相邻图像对进行评分的数据结构
  typedef hs::sfm::essential::EMatrix5PointsRansacRefiner<Scalar>
          EMatrixRansacRefiner;		//本质矩阵的RANSAC优化
  typedef typename EMatrixRansacRefiner::IndexSet IndexSet;	//来自 Ransac 模板类, 取其点类型的迭代器差距为Index, 再作为vector模板参数得到 IndexSet.
  typedef hs::sfm::essential::EMatrix5PointsCalculator<Scalar>
          EMatrixCalculator;		//本质矩阵计算
  typedef typename EMatrixCalculator::HKey HKey;
  typedef typename EMatrixCalculator::HKeyPair HKeyPair;
  typedef typename EMatrixCalculator::HKeyPairContainer HKeyPairContainer;
  typedef typename IntrinsicParams::KMatrix KMatrix;	//相机的K矩阵
  typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;
  typedef typename EMatrixCalculator::EMatrix EMatrix;	//本质矩阵类型定义
  typedef hs::sfm::essential::EMatrixExtrinsicParamsPointsCalculator<Scalar>
          ExtrinsicParamsPointsCalculator;

  struct ImagePairSize	//一个相邻图像对的key_pair尺寸容器, 以key_pair数量为键进行比较
  {
    size_t image0;
    size_t image1;
    size_t number_of_key_pairs;

    bool operator < (const ImagePairSize& other) const
    {
      return number_of_key_pairs < other.number_of_key_pairs;
    }
  };

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
                 hs::sfm::TrackContainer& tracks,
                 size_t& best_identity_id,
                 size_t& best_relative_id,
                 ExtrinsicParams& relative_extrinsic_params,
                 PointContainer& points,
                 hs::sfm::ObjectIndexMap& track_point_map,
                 hs::sfm::ViewInfoIndexer& view_info_indexer) const
  {
#if 1
    std::vector<ImagePairSize> image_pair_sizes;
    auto image_pair_itr = matches.begin();
    auto image_pair_itr_end = matches.end();
    for (; image_pair_itr != image_pair_itr_end; ++image_pair_itr)	//从match集合中创建ImagePairSize对象集合
    {
      ImagePairSize image_pair_size;
      image_pair_size.image0 = image_pair_itr->first.first;
      image_pair_size.image1 = image_pair_itr->first.second;
      image_pair_size.number_of_key_pairs = image_pair_itr->second.size();
      image_pair_sizes.push_back(image_pair_size);
    }
    std::sort(image_pair_sizes.begin(), image_pair_sizes.end());	//对ImagePairSize对象集合进行排序
    auto itr_image_pair_size = image_pair_sizes.rbegin();	//ImagePairSize对象集合开始逆序遍历
    auto itr_image_pair_size_end = image_pair_sizes.rend();
    best_identity_id = std::numeric_limits<size_t>::max();
    best_relative_id = std::numeric_limits<size_t>::max();
    IndexSet estimated_inlier_indices;	//Index的向量集合, 其中Index类型来自Ransac类, 即vector<Point>的difference_type类型(对于vector通常是size_t).
    PointContainer points_essential;
    for (; itr_image_pair_size != itr_image_pair_size_end;
         ++itr_image_pair_size)	//遍历相邻图像对
    {
      if (itr_image_pair_size->number_of_key_pairs <
          min_number_of_pair_matches_) break;	//判断是否小于最小匹配点要求, 若不满足则可跳出循环不再处理.
      hs::sfm::ImagePair image_pair(itr_image_pair_size->image0,
                                    itr_image_pair_size->image1);
      auto itr_key_pairs = matches.find(image_pair);
      auto key_pair_itr = itr_key_pairs->second.begin();
      auto key_pair_itr_end = itr_key_pairs->second.end();
      HKeyPairContainer key_pairs;
      size_t intrinsic_id_left = image_intrinsic_map[image_pair.first];		//根据图像ID获取其对应的内参数矩阵ID, 下同
      size_t intrinsic_id_right = image_intrinsic_map[image_pair.second];	
      const IntrinsicParams& intrinsic_params_left =
        intrinsic_params_set[intrinsic_id_left];	//获取内参数矩阵, 下同
      const IntrinsicParams& intrinsic_params_right =
        intrinsic_params_set[intrinsic_id_right];
      KMatrix K_left_inverse = intrinsic_params_left.GetKMatrix().inverse();	//获取相机K矩阵的逆
      KMatrix K_right_inverse = intrinsic_params_right.GetKMatrix().inverse();
      for (; key_pair_itr != key_pair_itr_end; ++key_pair_itr)
      {
        size_t key_left_id = key_pair_itr->first;
        size_t key_right_id = key_pair_itr->second;
        HKeyPair key_pair;	//HKeyPair 是EIGEN_VECTOR(Scalar, 3)的pair
        key_pair.first.segment(0, 2) =
          image_keysets[image_pair.first][key_left_id];		//从ImageKeys类型的对象中取点坐标, 类型为EIGEN_VECTOR(Scalar, 2), 下同.
        key_pair.first[2] = Scalar(1);	//齐次坐标, 第三位置1.
        key_pair.second.segment(0, 2) =
          image_keysets[image_pair.second][key_right_id];
        key_pair.second[2] = Scalar(1);
        key_pair.first = K_left_inverse * key_pair.first;		//取得图像点的归一化坐标, 下同.
        key_pair.second = K_right_inverse * key_pair.second;
        key_pairs.push_back(key_pair);
      }
      EMatrixRansacRefiner ransac_refiner;
      HKeyPairContainer key_pairs_refined;
      EMatrix e_matrix;
      estimated_inlier_indices.clear();
      if (ransac_refiner(key_pairs, 8 / intrinsic_params_left.focal_length(),
                         key_pairs_refined, estimated_inlier_indices,
                         e_matrix))	//若返回非0值,则表示RANSAC求解无效, 跳过 
      {
        continue;
      }

      //EMatrixCalculator ematrix_calculator;
      //if (ematrix_calculator(key_pairs_refined, e_matrix) != 0)
      //{
      //  continue;
      //}

      //通过E矩阵计算影像对的相对外方位元素
      ExtrinsicParamsPointsCalculator extrinsic_points_calculator;
      if (extrinsic_points_calculator(e_matrix,
                                      key_pairs_refined,
                                      relative_extrinsic_params,
                                      points_essential) != 0)	//在 求解成功 情况下返回0
      {
        continue;
      }

      best_identity_id = image_pair.first;	//找到最佳匹配image_pair, 准备输出.
      best_relative_id = image_pair.second;
      break;
    }
    if (best_identity_id != std::numeric_limits<size_t>::max() &&
        best_relative_id != std::numeric_limits<size_t>::max())		//若图像的ID有效, 则继续
    {
      auto itr_image_pair =
        matches.find(std::make_pair(best_identity_id, best_relative_id));	//取得对应的match
      if (itr_image_pair == matches.end())
      {
        return -1;
      }

      for (size_t i = 0; i < itr_image_pair->second.size(); i++)	//遍历该match的所有点对key_pair
      {
        size_t key_id_identity = itr_image_pair->second[i].first;
        size_t key_id_relative = itr_image_pair->second[i].second;
        ViewInfo* view_info_identity =
          view_info_indexer.GetViewInfoByImageKey(best_identity_id,
                                                     key_id_identity);	//构造key_pair对应的ViewInfo
        ViewInfo* view_info_relative =
          view_info_indexer.GetViewInfoByImageKey(best_relative_id,
                                                     key_id_relative);
        if (view_info_identity != nullptr)
        {
          view_info_identity->is_blunder = true;	//预设该成员的值?
        }
        if (view_info_relative != nullptr)
        {
          view_info_relative->is_blunder = true;
        }
      }

      points.clear();	//准备输出
      for (size_t i = 0; i < estimated_inlier_indices.size(); i++)	//遍历最佳匹配图像对的内点集合, 将其对应的match中, 关系正确的点输出.
      {
        size_t key_pair_id = size_t(estimated_inlier_indices[i]);
        size_t key_id_identity = itr_image_pair->second[key_pair_id].first;
        size_t key_id_relative = itr_image_pair->second[key_pair_id].second;
        ViewInfo* view_info_identity =
          view_info_indexer.GetViewInfoByImageKey(best_identity_id,
                                                     key_id_identity);
        ViewInfo* view_info_relative =
          view_info_indexer.GetViewInfoByImageKey(best_relative_id,
                                                     key_id_relative);
        if (view_info_identity != nullptr &&
            view_info_relative != nullptr)
        {
          view_info_identity->is_blunder = false;
          view_info_relative->is_blunder = false;
          size_t track_id = view_info_identity->track_id;	//记录track ID
          if (track_id == view_info_relative->track_id)		//记录并输出.
          {
            track_point_map[track_id] = points.size();
            points.push_back(points_essential[i]);
          }
        }
      }
//#if 1
//      size_t number_of_points = 0;
//#endif
//      for (size_t i = 0; i < tracks.size(); i++)
//      {
//        size_t number_of_views = tracks[i].size();
//        size_t key_id_identity = std::numeric_limits<size_t>::max();
//        size_t key_id_relative = std::numeric_limits<size_t>::max();
//        for (size_t j = 0; j < number_of_views; j++)
//        {
//          size_t image_id = tracks[i][j].first;
//          size_t key_id = tracks[i][j].second;
//          if (image_id == best_identity_id)
//          {
//            key_id_identity = key_id;
//          }
//          if (image_id == best_relative_id)
//          {
//            key_id_relative = key_id;
//          }
//        }
//        if (key_id_identity != std::numeric_limits<size_t>::max() &&
//            key_id_relative != std::numeric_limits<size_t>::max())
//        {
//          auto itr_key_pair =
//            key_pair_indexer.find(std::make_pair(key_id_identity,
//                                                 key_id_relative));
//          if (itr_key_pair != key_pair_indexer.end())
//          {
//            track_point_map[i] = itr_key_pair->second;
//#if 1
//            std::cout << "key_id_identity:" << key_id_identity << "\n";
//            std::cout << "key_id_relative:" << key_id_relative << "\n";
//            std::cout << "point_id:" << itr_key_pair->second << "\n";
//            number_of_points++;
//#endif
//          }
//          else
//          {
//            for (size_t j = 0; j < number_of_views; j++)
//            {
//              ViewInfo& view_info = view_info_indexer.GetViewInfoByTrackImage(
//                                      i, tracks[i][j].first);
//              view_info.is_blunder = true;
//            }
//          }
//        }
//      }

//#if 1
//      std::cout << "number_of_points:" << number_of_points << "\n";
//      std::cout << "points.size():" << points.size() << "\n";
//#endif

      return 0;
    }
    else
    {
      return -1;
    }
#else
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

    return 0;
#endif
  }

private:
  size_t min_number_of_pair_matches_;
};

}
}
}

#endif
