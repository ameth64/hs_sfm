#ifndef _HS_SFM_SFM_PIPELINE_IMAGE_EXPANDOR_HPP_
#define _HS_SFM_SFM_PIPELINE_IMAGE_EXPANDOR_HPP_

#include <map>
#include <algorithm>
#include <cmath>

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/projective_functions.hpp"
#include "hs_sfm/projective/single_camera_params_maximum_likelihood_estimator.hpp"
#include "hs_sfm/projective/pmatrix_dlt_ransac_refiner.hpp"

#define DEBUG_TMP 1
#if DEBUG_TMP
#include <iostream>
#endif

namespace hs
{
namespace sfm
{
namespace pipeline
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
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

  typedef ObjectIndexMap TrackPointMap;
  typedef ObjectIndexMap ImageIntrinsicMap;
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
                  const ImageIntrinsicMap& image_intrinsic_map,
                  const PointContainer& points,
                  const TrackPointMap& track_point_map,
                  const ImageViewTracksContainer& image_view_tracks_set,
                  const ImageExtrinsicMap& image_extrinsic_map,
                  ViewInfoIndexer& view_info_indexer,
                  ExtrinsicParamsContainer& new_extrinsic_params_set,
                  std::vector<size_t>& new_image_ids) const
  {
    size_t number_of_images = image_keysets.size();
    if (number_of_images != image_intrinsic_map.Size())
    {
      return -1;
    }

    //查找拍摄到当前点云的其他影像，设匹配数最多的影像为n，
    //选取所有匹配数大于0.75n的影像，
    //若匹配数最多的影像的匹配数没超过阈值，则函数执行不成功。

    size_t max_matches_image_id = 0;
    size_t max_number_of_view_tracks = 0;
    std::map<size_t, size_t> number_of_view_tracks_map;
    std::map<size_t, ImageViewTracks> image_view_tracks_map;
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

      number_of_view_tracks_map[i] = number_of_view_tracks;
      image_view_tracks_map[i] = image_view_tracks;

      if (number_of_view_tracks > max_number_of_view_tracks)
      {
        max_number_of_view_tracks = number_of_view_tracks;
        max_matches_image_id = i;
      }
    }
    if (max_number_of_view_tracks < add_new_image_matches_threshold_)
    {
      return -1;
    }

    auto itr_number_of_view_tracks = number_of_view_tracks_map.begin();
    auto itr_number_of_view_tracks_end = number_of_view_tracks_map.end();
    auto itr_image_view_tracks = image_view_tracks_map.begin();
    auto itr_image_view_tracks_end = image_view_tracks_map.end();
    for (; itr_number_of_view_tracks != itr_number_of_view_tracks_end &&
           itr_image_view_tracks != itr_image_view_tracks_end;
         ++itr_number_of_view_tracks, ++itr_image_view_tracks)
    {
      size_t number_of_view_tracks = itr_number_of_view_tracks->second;
      if (number_of_view_tracks < size_t(0.75 * max_number_of_view_tracks))
      {
        continue;
      }
      const ImageViewTracks& image_view_tracks = itr_image_view_tracks->second;
      size_t image_id = itr_number_of_view_tracks->first;
      //获取三维点与二维点的对应
      CorrespondenceContainer coarse_correspondences;
      ObjectIndexMap key_track_map(number_of_view_tracks);
      for (size_t i = 0; i < number_of_view_tracks; i++)
      {
        size_t track_id = image_view_tracks[i];
        if (!track_point_map.IsValid(track_id))
        {
          return -1;
        }
        size_t point_id = track_point_map[track_id];

        ViewInfo& view_info =
          view_info_indexer.GetViewInfoByTrackImage(track_id, image_id);
        view_info.is_blunder = true;

        Correspondence correspondence;
        correspondence.first = image_keysets[image_id][view_info.key_id];
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
        continue;
      }

      size_t number_of_inliers = inlier_indices.size();
      if (number_of_inliers < add_new_image_matches_threshold_)
      {
        continue;
      }

      //使用剔除粗差的点计算外方位元素
      CameraParamsMLEstimator camera_ml_estimator;
      KeyCovariance key_covariance = KeyCovariance::Identity();
      IntrinsicParams intrinsic_params_estimate;
      size_t intrinsic_id = image_intrinsic_map[image_id];
      ExtrinsicParams new_extrinsic_params;
      IndexSet mle_inlier_indices = inlier_indices;
      CorrespondenceContainer mle_correspondences = refined_correspondences;
      bool estimate_success = false;
      while (1)
      {
        if (camera_ml_estimator(mle_correspondences,
                                key_covariance,
                                intrinsic_params_set[intrinsic_id],
                                Scalar(1),
                                Scalar(0.01),
                                Scalar(1),
                                Scalar(1),
                                Scalar(1e-4),
                                Scalar(1e-5),
                                Scalar(1e-5),
                                Scalar(1e-5),
                                Scalar(1e-7),
                                Scalar(1e-7),
                                intrinsic_params_estimate,
                                new_extrinsic_params) == 0)
        {
          estimate_success = true;
          std::vector<Scalar> errors(mle_correspondences.size());
          for (size_t i = 0; i < mle_correspondences.size(); i++)
          {
            Key key =
              hs::sfm::ProjectiveFunctions<Scalar>::WorldPointProjectToImageKey(
                intrinsic_params_estimate, new_extrinsic_params,
                mle_correspondences[i].second);

            errors[i] = (key - mle_correspondences[i].first).norm();
          }

          size_t nth = size_t(0.95 * double(mle_correspondences.size()));
          std::vector<Scalar> errors_copy = errors;
          std::nth_element(errors_copy.begin(),
                           errors_copy.begin() + nth,
                           errors_copy.end());
          Scalar nth_median = errors_copy[nth];
          Scalar threshold = Scalar(2.4) * nth_median;
#if DEBUG_TMP
          std::cout<<"threshold:"<<threshold<<"\n";
#endif
          threshold = std::max(std::min(threshold, Scalar(24)), Scalar(12));
          CorrespondenceContainer mle_correspondences_next;
          IndexSet mle_inlier_indices_next;
          for (size_t i = 0; i < errors.size(); i++)
          {
            if (errors[i] < threshold)
            {
              mle_correspondences_next.push_back(mle_correspondences[i]);
              mle_inlier_indices_next.push_back(mle_inlier_indices[i]);
            }
#if DEBUG_TMP
            else
            {
              std::cout<<"Add image error:"<<errors[i]
                       <<" threshold:"<<threshold<<"\n";
            }
#endif
          }

          if (mle_correspondences_next.size() == mle_correspondences.size())
          {
            break;
          }
          else if (mle_correspondences_next.empty())
          {
            estimate_success = false;
            break;
          }
          else
          {
            mle_correspondences.swap(mle_correspondences_next);
            mle_inlier_indices.swap(mle_inlier_indices_next);
          }
        }
        else
        {
          break;
        }
      }

      if (estimate_success)
      {
        for (size_t i = 0; i < mle_inlier_indices.size(); i++)
        {
          size_t track_id =
            key_track_map[mle_inlier_indices[i]];
          ViewInfo& view_info =
            view_info_indexer.GetViewInfoByTrackImage(track_id, image_id);
          view_info.is_blunder = false;
        }
        new_extrinsic_params_set.push_back(new_extrinsic_params);
        new_image_ids.push_back(image_id);
      }
    }

    if (new_extrinsic_params_set.empty())
    {
      return -1;
    }
    else
    {
      return 0;
    }
  }
private:
  size_t add_new_image_matches_threshold_;
  Scalar pmatrix_ransac_threshold_;
};

}
}
}

#endif
