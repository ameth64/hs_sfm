#ifndef _HS_SFM_INCREMENTAL_SCENE_SPANNER_HPP_
#define _HS_SFM_INCREMENTAL_SCENE_SPANNER_HPP_

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
class SceneSpanner
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
  typedef std::map<size_t, size_t> TrackPointMap;
  typedef std::map<size_t, size_t> ImageExtrinsicMap;

private:
  struct ViewInfo
  {
    size_t image_id;
    size_t key_id;
    bool is_outlier;
  };
  typedef std::vector<ViewInfo> RichTrack;
  typedef std::vector<RichTrack> RichTrackContainer;
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
  SceneSpanner(
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
                  TrackPointMap& track_point_map) const
  {
    size_t number_of_images = image_keysets.size();
    if (number_of_images != intrinsic_params_set.size())
    {
      return -1;
    }

    //构造RichTrackContainer
    //计算每个内参数对应的影像拍到的track
    ImageViewTracksContainer image_view_tracks_set(number_of_images);
    RichTrackContainer rich_tracks;
    size_t number_of_tracks = tracks.size();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();
      RichTrack rich_track;
      for (size_t j = 0; j < number_of_views; j++)
      {
        image_view_tracks_set[tracks[i][j].first].push_back(i);
        ViewInfo view;
        view.image_id = tracks[i][j].first;
        view.key_id = tracks[i][j].second;
        view.is_outlier = false;
        rich_track.push_back(view);
      }
      rich_tracks.push_back(rich_track);
    }

    ExtrinsicParams new_extrinsic_params;
    size_t new_image_id;
    AddNewImage(image_keysets,
                intrinsic_params_set,
                points,
                track_point_map,
                image_view_tracks_set,
                image_extrinsic_map,
                rich_tracks,
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
                  RichTrackContainer& rich_tracks,
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
      if (image_extrinsic_map.find(i) != image_extrinsic_map.end())
      {
        continue;
      }

      size_t number_of_view_tracks = 0;
      ImageViewTracks image_view_tracks;
      for (size_t j = 0; j < image_view_tracks_set[i].size(); j++)
      {
        if (track_point_map.find(image_view_tracks_set[i][j]) ==
            track_point_map.end())
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
    std::vector<std::pair<size_t, size_t> > key_track_map;
    for (size i = 0; i < max_number_of_view_tracks; i++)
    {
      size_t track_id = image_view_tracks_set[i];
      const RichTrack& rich_track = rich_tracks[track_id];
      auto itr_track_point = track_point_map.find(image_view_tracks_set[i]);
      if (itr_track_point == track_point_map.end())
      {
        return -1;
      }
      size_t point_id = itr_track_point->second;
      size_t number_of_views = rich_track.size();
      size_t key_id = std::numeric_limits<size_t>::max();
      size_t view_id = 0;
      for (size_t j = 0; j < number_of_views; j++)
      {
        if (rich_track[j].image_id == new_image_id)
        {
          key_id = rich_track[j].key_id;
          view_id = j;
        }
      }
      if (key_id == std::numeric_limits<size_t>::max())
      {
        return -1;
      }

      Correspondence correspondence;
      correspondence.first = image_keysets[new_image_id][key_id];
      correspondence.second = points[point_id];
      coarse_correspondences.push_back(correspondence);
      key_track_map.push_back(std::pair<size_t, size_t>(track_id, view_id));
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
      
    }

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