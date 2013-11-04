#ifndef _HS_SFM_INCREMENTAL_SCENE_EXPANDOR_HPP_
#define _HS_SFM_INCREMENTAL_SCENE_EXPANDOR_HPP_

#include <vector>
#include <map>

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/projective/pmatrix_dlt_calculator.hpp"
#include "hs_sfm/projective/pmatrix_dlt_ransac_refiner.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar>
class SceneExpandor
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef ImageKeys<Scalar> ImageKeyset;
  typedef EIGEN_STD_VECTOR(ImageKeyset) ImageKeysetContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;
  typedef std::vector<size_t> TrackPointMap;
  typedef std::vector<size_t> ImageExtrinsicMap;

  class ObjectIndexMap
  {
  public:
    static const size_t invalid_value_ = std::numeric_limits<size_t>::max();
    ObjectIndexMap(){}
    ObjectIndexMap(size_t number_of_objects)
      : mapper_(number_of_objects, invalid_value) {}
  public:
    size_t GetMappedId(size_t object_id) const
    {
      return mapper_[object_id];
    }

    void SetObjectId(size_t object_id, size_t mapped_id)
    {
      mapper_[object_id] = mapped_id;
    }

    size_t operator[] (size_t object_id) const
    {
      return mapper_[object_id];
    }

    size_t& operator[] (size_t object_id)
    {
      return mapper_[object_id];
    }

    bool IsValid(size_t object_id) const
    {
      return mapper_[object_id] != invalid_value_;
    }

    void Resize(size_t number_of_objects)
    {
      mapper_.resize(number_of_objects, invalid_value_);
    }
  private:
    std::vector<size_t> mapper_;
  };

  typedef ObjectIndexMap TrackPointMap;
  typedef ObjectIndexMap ImageExtrinsicMap;

private:
  struct ViewInfo
  {
    size_t track_id;
    size_t image_id;
    size_t key_id;
    bool is_blunder;
  };

  class ViewInfoIndexer
  {
  public:
    void SetViewInfoByTracks(const hs::sfm::TrackContainer& tracks,
                             const TrackPointMap& track_point_map,
                             const ImageExtrinsicMap& image_extrinsic_map)
    {
      views_info_.cleaa();
      size_t number_of_tracks = tracks.size();
      for (size_t i = 0; i < number_of_tracks; i++)
      {
        size_t number_of_views = tracks[i].size();
        for (size_t j = 0; j < number_of_views; j++)
        {
          ViewInfo view_info;
          view_info.track_id = i;
          view_info.image_id = tracks[i][j].first;
          view_info.key_id = tracks[i][j].second;
          view_info.is_blunder = (track_point_map.IsValid(i) &&
                                  image_extrinsic_map[view_info.image_id]);
          views_info.push_back(view_info);
        }
      }
    }

    const ViewInfo& GetViewInfoByTrackImage(size_t track_id,
                                            size_t image_id) const
    {
      return views_info[track_image_index[std::make_pair(track_id, image_id)]];
    }

    ViewInfo& GetViewInfoByTrackImage(size_t track_id,
                                      size_t image_id)
    {
      return views_info[track_image_index[std::make_pair(track_id, image_id)]];
    }

    const ViewInfo& GetViewInfoByImageKey(size_t image_id,
                                        size_t key_id) const
    {
      return views_info[image_key_index[std::make_pair(image_id, key_id)]];
    }

    ViewInfo& GetViewInfoByImageKey(size_t image_id,
                                  size_t key_id)
    {
      return views_info[image_key_index[std::make_pair(image_id, key_id)]];
    }

  private:
    std::vector<ViewInfo> views_info;
    std::map<std::pair<size_t, size_t>, size_t> track_image_index;
    std::map<std::pair<size_t, size_t>, size_t> image_key_index;
  };

  typedef std::vector<size_t> ImageViewTracks;
  typedef std::vector<ImageViewTracks> ImageViewTracksContainer;
  typedef hs::sfm::projective::PMatrixDLTCalculator<Scalar>
          PMatrixCalculator;
  typedef hs::sfm::projective::PMatrixDLTRansacRefiner<Scalar>
          PMatrixRasacRefiner;
  typedef typename PMatrixCalculator::Key Key;
  typedef typename PMatrixCalculator::Correspondence Correspondence;
  typedef typename PMatrixCalculator::CorrespondenceContainer
                   CorrespondenceContainer;
  typedef typename PMatrixCalculator::PMatrix PMatrix;
  typedef typename PMatrixRasacRefiner::IndexSet IndexSet;

public:
  SceneExpandor(
    size_t add_new_image_matches_threshold,
    Scalar pmatrix_ransac_threshold)
    : add_new_image_matches_threshold_(add_new_image_matches_threshold),
      pmatrix_ransac_threshold_(pmatrix_ransac_threshold){}

  Err operator() (const ImageKeysetContainer& image_keysets,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const TrackContainer& tracks,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  ImageExtrinsicMap& image_extrinsic_map,
                  PointContainer& points,
                  TrackPointMap& track_point_map,
                  ViewInfoIndexer& view_info_indexer) const
  {
    size_t number_of_images = image_keysets.size();
    if (number_of_images != intrinsic_params_set.size())
    {
      return -1;
    }

    //构造view info indexer
    //计算每个内参数对应的影像拍到的track
    ImageViewTracksContainer image_view_tracks_set(number_of_images);
    size_t number_of_tracks = tracks.size();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();
      for (size_t j = 0; j < number_of_views; j++)
      {
        image_view_tracks_set[tracks[i][j].first].push_back(i);
      }
    }
    view_info_indexer.SetViewInfoByTracks(tracks,
                                          track_point_map,
                                          image_extrinsic_map);

    ExtrinsicParams new_extrinsic_params;
    size_t new_image_id;
    AddNewImage(image_keysets,
                intrinsic_params_set,
                points,
                track_point_map,
                image_view_tracks_set,
                image_extrinsic_map,
                view_info_indexer,
                new_extrinsic_params,
                new_image_id);

    return 0;
  }

private:
  Err AddNewImage(const ImageKeysetContainer& image_keysets,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const PointContainer& points,
                  const TrackPointMap& track_point_map,
                  const ImageViewTracksContainer& image_view_tracks_set,
                  const ImageExtrinsicMap& image_extrinsic_map,
                  ViewInfoIndexer& view_info_indexer,
                  ExtrinsicParams& new_extrinsic_params,
                  size_t& new_image_id) const
  {
    size_t number_of_images = image_keysets.size();
    if (number_of_images != intrinsic_params_set.size())
    {
      return -1;
    }

    //查找拍摄到当前点云的其他影像，选取匹配数最多的影像，
    //若匹配数最多的影像的匹配数没超过阈值，则函数执行不成功。
    size_t max_matches_image_id = 0;
    size_t max_number_of_view_tracks = 0;
    ImageViewTracks max_image_view_tracks;
    for (size_t i = 0; i < number_of_images; i++)
    {
      //已加入的影像排除在外
      if (image_extrinsic_map.IsValid(i))
      {
        continue;
      }

      size_t number_of_view_tracks = 0;
      ImageViewTracks image_view_tracks;
      for (size_t j = 0; j < image_view_tracks_set[i].size(); j++)
      {
        if (!track_point_map.IsValid(image_view_tracks_set[i][j]))
        {
          continue;
        }

        image_view_tracks.push_back(image_view_tracks_set[i][j]);
        number_of_view_tracks++;
      }

      if (number_of_view_tracks > max_number_of_view_tracks)
      {
        max_number_of_view_tracks = number_of_view_tracks;
        max_matches_image_id = i;
        max_image_view_tracks.swap(image_view_tracks);
      }
    }

    if (max_number_of_view_tracks < add_new_image_matches_threshold_)
    {
      return -1;
    }

    new_image_id = max_matches_image_id;

    //获取三维点与二维点的对应
    CorrespondenceContainer coarse_correspondences;
    ObjectIndexMap key_track_map(max_number_of_view_tracks);
    for (size_t i = 0; i < max_number_of_view_tracks; i++)
    {
      size_t track_id = image_view_tracks_set[new_image_id][i];
      if (!track_point_map.IsValid(track_id))
      {
        return -1;
      }

      const ViewInfo view_info =
        view_info_indexer.GetViewInfoByTrackImage(track_id, new_image_id);

      Correspondence correspondence;
      correspondence.first = image_keysets[new_image_id][view_info.key_id];
      correspondence.second = points[point_id];
      coarse_correspondences.push_back(correspondence);
      key_track_map[i] = track_id;
    }

    //Ransac剔除错误的对应
    IndexSet inlier_indices;
    PMatrixRasacRefiner ransac_refiner;
    CorrespondenceContainer refined_correspondences;
    if (ransac_refiner(coarse_correspondences,
                       pmatrix_ransac_threshold_,
                       refined_correspondences,
                       inlier_indices) != 0)
    {
      return -1;
    }
    size_t number_of_inliers = inlier_indices.size();
    for (size_t i = 0; i < number_of_inliers; i++)
    {
      size_t track_id = key_track_map[inlier_indices[i]];
      ViewInfo& view_info =
        view_info_indexer.GetViewInfoByTrackImage(track_id, new_image_id);
      view_info.is_blunder = false;
    }

    //使用剔除粗差的点计算P矩阵
    PMatrixDLTCalculator calculator;
    PMatrix p_matrix;
    if (calculator(refined_correspondences, p_matrix) != 0) return -1;

    //由P矩阵获取影像外方位元素

    return 0;
  }

  Err TriangulateNewPoints(const ImageKeysetContainer& image_keysets,
                           const IntrinsicParamsContainer& intrinsic_params_set,
                           const TrackContainer& tracks,
                           const ExtrinsicParams& extrinsic_params_set,
                           const ImageExtrinsicMap& image_extrinsic_map,
                           PointContainer& points,
                           TrackPointMap& track_point_map) const
  {
    return 0;
  }

  Err BundleAdjustment(const ImageKeysetContainer& image_keysets,
                       const IntrinsicParamsContainer& intrinsic_params_set,
                       const TrackContainer& tracks,
                       const ImageExtrinsicMap& image_extrinsic_map,
                       ExtrinsicParamsContainer& extrinsic_params_set,
                       PointContainer& points,
                       TrackPointMap& point_track_map) const
  {
    return 0;
  }

private:
  size_t add_new_image_matches_threshold_;
  Scalar pmatrix_ransac_threshold_;
};

}
}
}

#endif
