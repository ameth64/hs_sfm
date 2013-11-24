#ifndef _HS_SFM_INCREMENTAL_INCREMENTAL_HPP_
#define _HS_SFM_INCREMENTAL_INCREMENTAL_HPP_

#include <map>
#include <limits>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"
#include "hs_sfm/incremental/best_pair_selector.hpp"
#include "hs_sfm/incremental/scene_expandor.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar>
class IncrementalSFM
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

private:
  typedef BestPairSelector<Scalar> BestPairSelector;
  typedef typename ExtrinsicParams::Rotation Rotation;

public:
  IncrementalSFM()
    : min_number_of_pair_matches_(100),
      add_new_image_matches_threshold_(8),
      min_triangulate_views_(2) {}

  IncrementalSFM(size_t min_number_of_pair_matches,
                 size_t add_new_image_matches_threshold,
                 size_t min_triangulate_views)
    : min_number_of_pair_matches_(min_number_of_pair_matches),
      add_new_image_matches_threshold_(add_new_image_matches_threshold),
      min_triangulate_views_(min_triangulate_views) {}

  Err operator() (const IntrinsicParamsContainer& intrinsic_params_set,
                  const hs::sfm::MatchContainer& matches,
                  const KeysetContainer& keysets,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  hs::sfm::ObjectIndexMap& image_extrinsic_map,
                  PointContainer& points,
                  hs::sfm::ObjectIndexMap& track_point_map,
                  hs::sfm::ViewInfoIndexer& view_info_indexer) const
  {
    BestPairSelector selector(min_number_of_pair_matches_);
    size_t best_identity_id;
    size_t best_relative_id;
    ExtrinsicParams extrinsic_params_relative;
    PointContainer points_best_pair;
    if (selector(keysets,
                 matches,
                 intrinsic_params_set,
                 best_identity_id,
                 best_relative_id,
                 extrinsic_params_relative,
                 points_best_pair) != 0) return -1;

    //旋转向量为0时（即没有旋转），bundle adjustment时无法计算jacobian矩阵。
    //因此需改变初始的相机朝向。
    Rotation rotation_extra;
    rotation_extra[0] = Scalar(0);
    rotation_extra[1] = Scalar(0);
    rotation_extra[2] = Scalar(3.141592653) / 180;

    extrinsic_params_relative.rotation() =
      extrinsic_params_relative.rotation() * rotation_extra.Inverse();

    extrinsic_params_relative.position() =
      rotation_extra * extrinsic_params_relative.position();

    ExtrinsicParams extrinsic_params_identity;
    extrinsic_params_identity.rotation() = rotation_extra.Inverse();
    extrinsic_params_identity.position().setZero();

    size_t number_of_points = points_best_pair.size();
    for (size_t i = 0; i < number_of_points; i++)
    {
      Point point = points_best_pair[i];
      points_best_pair[i] = rotation_extra * point;
    }

    size_t number_of_images = intrinsic_params_set.size();
    extrinsic_params_set.clear();
    extrinsic_params_set.push_back(extrinsic_params_identity);
    extrinsic_params_set.push_back(extrinsic_params_relative);
    image_extrinsic_map.Resize(number_of_images);
    image_extrinsic_map[best_identity_id] = 0;
    image_extrinsic_map[best_relative_id] = 1;

    TrackContainer tracks;
    hs::sfm::MatchesTracksConvertor matches_track_convertor;
    if (matches_track_convertor(matches, tracks) != 0)
    {
      return -1;
    }
    size_t number_of_tracks = tracks.size();
    track_point_map.Resize(number_of_tracks);
    std::map<std::pair<size_t, size_t>, size_t> key_pair_indexer;
    auto itr_image_pair = matches.find(std::make_pair(best_identity_id,
                                                      best_relative_id));
    if (itr_image_pair == matches.end())
    {
      return -1;
    }
    size_t number_of_key_pairs = itr_image_pair->second.size();
    for (size_t i = 0; i < number_of_key_pairs; i++)
    {
      key_pair_indexer[std::make_pair(itr_image_pair->second[i].first,
                                      itr_image_pair->second[i].second)] = i;
    }
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();
      size_t key_id_identity = std::numeric_limits<size_t>::max();
      size_t key_id_relative = std::numeric_limits<size_t>::max();
      for (size_t j = 0; j < number_of_views; j++)
      {
        size_t image_id = tracks[i][j].first;
        size_t key_id = tracks[i][j].second;
        if (image_id == best_identity_id)
        {
          key_id_identity = key_id;
        }
        if (image_id == best_relative_id)
        {
          key_id_relative = key_id;
        }
      }
      if (key_id_identity != std::numeric_limits<size_t>::max() &&
          key_id_relative != std::numeric_limits<size_t>::max())
      {
        auto itr_key_pair =
          key_pair_indexer.find(std::make_pair(key_id_identity,
                                               key_id_relative));
        if (itr_key_pair == key_pair_indexer.end())
        {
          return -1;
        }
        track_point_map[i] = itr_key_pair->second;
      }
    }

    points.swap(points_best_pair);
    SceneExpandor<Scalar> expandor(add_new_image_matches_threshold_,
                                   Scalar(4),
                                   min_triangulate_views_,
                                   Scalar(4));
    if (expandor(keysets,
                 intrinsic_params_set,
                 tracks,
                 extrinsic_params_set,
                 image_extrinsic_map,
                 points,
                 track_point_map,
                 view_info_indexer) != 0)
    {
      return -1;
    }

    return 0;
  }

private:
  size_t min_number_of_pair_matches_;
  size_t add_new_image_matches_threshold_;
  size_t min_triangulate_views_;
};

}
}
}

#endif
