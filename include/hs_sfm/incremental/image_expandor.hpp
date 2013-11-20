﻿#ifndef _HS_SFM_INCREMENTAL_IMAGE_EXPANDOR_HPP_
#define _HS_SFM_INCREMENTAL_IMAGE_EXPANDOR_HPP_

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/projective/single_camera_params_maximum_likelihood_estimator.hpp"
#include "hs_sfm/projective/pmatrix_dlt_ransac_refiner.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar>
class ImageExpandor
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

  typedef ObjectIndexMap TrackPointMap;
  typedef ObjectIndexMap ImageExtrinsicMap;
  typedef std::vector<size_t> ImageViewTracks;
  typedef std::vector<ImageViewTracks> ImageViewTracksContainer;

private:
  typedef
    hs::sfm::projective::SingleCameraParamsMaximumLikelihoodEstimator<Scalar>
      CameraParamsMLEstimator;
  typedef typename CameraParamsMLEstimator::Correspondence Correspondence;
  typedef typename CameraParamsMLEstimator::CorrespondenceContainer
                   CorrespondenceContainer;
  typedef typename CameraParamsMLEstimator::KeyCovariance KeyCovariance;

  typedef hs::sfm::projective::PMatrixDLTRansacRefiner<Scalar>
          PMatrixRasacRefiner;
  typedef typename PMatrixRasacRefiner::IndexSet IndexSet;

public:
  ImageExpandor()
    : add_new_image_matches_threshold_(8),
      pmatrix_ransac_threshold_(4.0) {}

  ImageExpandor(size_t add_new_image_matches_threshold,
                Scalar pmatrix_ransac_threshold)
    : add_new_image_matches_threshold_(add_new_image_matches_threshold),
      pmatrix_ransac_threshold_(pmatrix_ransac_threshold) {}

  Err operator() (const ImageKeysetContainer& image_keysets,
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
      size_t track_id = max_image_view_tracks[i];
      if (!track_point_map.IsValid(track_id))
      {
        return -1;
      }
      size_t point_id = track_point_map[track_id];

      ViewInfo& view_info =
        view_info_indexer.GetViewInfoByTrackImage(track_id, new_image_id);
      view_info.is_blunder = true;

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

    //使用剔除粗差的点计算外方位元素
    CameraParamsMLEstimator camera_ml_estimator;
    KeyCovariance key_covariance = KeyCovariance::Identity();
    IntrinsicParams intrinsic_params_estimate;
    if (camera_ml_estimator(refined_correspondences,
                            key_covariance,
                            intrinsic_params_set[new_image_id],
                            Scalar(1),
                            Scalar(0.01),
                            Scalar(1),
                            Scalar(1),
                            Scalar(1e-4),
                            intrinsic_params_estimate,
                            new_extrinsic_params) != 0)
    {
      return -1;
    }

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
