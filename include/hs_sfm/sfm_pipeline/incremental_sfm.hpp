#ifndef _HS_SFM_SFM_PIPELINE_INCREMENTAL_HPP_
#define _HS_SFM_SFM_PIPELINE_INCREMENTAL_HPP_

#include <map>
#include <limits>
#include <algorithm>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_progress/progress_utility/progress_manager.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"
#include "hs_sfm/sfm_pipeline/best_pair_selector.hpp"
#include "hs_sfm/sfm_pipeline/scene_expandor.hpp"

#define DEBUG_TMP 1
#if DEBUG_TMP
#include "hs_sfm/sfm_utility/debug_tmp.hpp"
#endif

namespace hs
{
namespace sfm
{
namespace pipeline
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

#if DEBUG_TMP
  typedef DebugTrue<Scalar> DebugTrueType;
#endif

private:
  typedef BestPairSelector<Scalar> BestPairSelector;
  typedef typename ExtrinsicParams::Rotation Rotation;

public:
  IncrementalSFM()
    : min_number_of_pair_matches_(100),
      add_new_image_matches_threshold_(8),
      min_triangulate_views_(2),
      number_of_threads_(1) {}

  IncrementalSFM(size_t min_number_of_pair_matches,
                 size_t add_new_image_matches_threshold,
                 size_t min_triangulate_views,
                 size_t number_of_threads)
    : min_number_of_pair_matches_(min_number_of_pair_matches),
      add_new_image_matches_threshold_(add_new_image_matches_threshold),
      min_triangulate_views_(min_triangulate_views),
      number_of_threads_(number_of_threads) {}

  Err operator() (const ObjectIndexMap& image_intrinsic_map,
                  const hs::sfm::MatchContainer& matches,
                  const KeysetContainer& keysets,
                  IntrinsicParamsContainer& intrinsic_params_set,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  hs::sfm::ObjectIndexMap& image_extrinsic_map,
                  PointContainer& points,
                  TrackContainer& tracks,
                  hs::sfm::ObjectIndexMap& track_point_map,
                  hs::sfm::ViewInfoIndexer& view_info_indexer,
                  hs::progress::ProgressManager* progress_manager = NULL
#if DEBUG_TMP
                  , const DebugTrueType& debug_true = DebugTrueType()
#endif
                  ) const
  {
    hs::sfm::MatchesTracksConvertor matches_track_convertor;
    if (matches_track_convertor(matches, tracks) != 0)
    {
      return -1;
    }
    auto itr_track = tracks.begin();
    auto itr_track_end = tracks.end();
    for (; itr_track != itr_track_end; ++itr_track)
    {
      std::sort(itr_track->begin(), itr_track->end());
    }
    std::sort(tracks.begin(), tracks.end());
    MatchContainer matches_filtered;
    if (matches_track_convertor(tracks, matches_filtered) != 0)
    {
      return -1;
    }

    size_t number_of_tracks = tracks.size();
    track_point_map.Resize(number_of_tracks);

    //构造view info indexer
    view_info_indexer.SetViewInfoByTracks(tracks);

    BestPairSelector selector(min_number_of_pair_matches_);
    size_t best_identity_id;
    size_t best_relative_id;
    ExtrinsicParams extrinsic_params_relative;
    PointContainer points_best_pair;
    if (selector(keysets,
                 matches_filtered,
                 intrinsic_params_set,
                 image_intrinsic_map,
                 tracks,
                 best_identity_id,
                 best_relative_id,
                 extrinsic_params_relative,
                 points_best_pair,
                 track_point_map,
                 view_info_indexer) != 0) return -1;

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
    size_t number_of_images = keysets.size();
    extrinsic_params_set.clear();
    extrinsic_params_set.push_back(extrinsic_params_identity);
    extrinsic_params_set.push_back(extrinsic_params_relative);
    image_extrinsic_map.Resize(number_of_images);
    image_extrinsic_map[best_identity_id] = 0;
    image_extrinsic_map[best_relative_id] = 1;

    points.swap(points_best_pair);
    SceneExpandor<Scalar> expandor(add_new_image_matches_threshold_,
                                   Scalar(16),
                                   min_triangulate_views_,
                                   Scalar(16),
                                   number_of_threads_);
    if (expandor(keysets,
                 image_intrinsic_map,
                 tracks,
                 intrinsic_params_set,
                 extrinsic_params_set,
                 image_extrinsic_map,
                 points,
                 track_point_map,
                 view_info_indexer,
                 progress_manager
#if DEBUG_TMP
                 , debug_true
#endif
                 ) != 0)
    {
      return -1;
    }

    return 0;
  }

private:
  size_t min_number_of_pair_matches_;
  size_t add_new_image_matches_threshold_;
  size_t min_triangulate_views_;
  size_t number_of_threads_;
};

}
}
}

#endif
