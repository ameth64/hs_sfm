﻿#ifndef _HS_SFM_SFM_PIPELINE_SCENE_EXPANDOR_HPP_
#define _HS_SFM_SFM_PIPELINE_SCENE_EXPANDOR_HPP_

#include <vector>
#include <map>

#include "hs_progress/progress_utility/progress_manager.hpp"

#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_pipeline/image_expandor.hpp"
#include "hs_sfm/sfm_pipeline/point_expandor.hpp"
#include "hs_sfm/sfm_pipeline/bundle_adjustment_optimizor.hpp"

namespace hs
{
namespace sfm
{
namespace pipeline
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

  typedef ObjectIndexMap TrackPointMap;
  typedef ObjectIndexMap ImageIntrinsicMap;
  typedef ObjectIndexMap ImageExtrinsicMap;

protected:
  typedef hs::sfm::pipeline::ImageExpandor<Scalar> ImageExpandor;
  typedef typename ImageExpandor::ImageViewTracks ImageViewTracks;
  typedef typename ImageExpandor::ImageViewTracksContainer
                   ImageViewTracksContainer;

  typedef hs::sfm::pipeline::PointExpandor<Scalar> PointExpandor;

  typedef hs::sfm::pipeline::BundleAdjustmentOptimizor<Scalar>
          BundleAdjustmentOptimizor;

public:
  SceneExpandor(
    size_t add_new_image_matches_threshold = 8,
    Scalar pmatrix_ransac_threshold = 4.0,
    size_t min_triangulate_views = 2,
    Scalar triangulate_error_threshold = 4.0,
    size_t number_of_threads = 1)
    : image_expandor_(add_new_image_matches_threshold,
                      pmatrix_ransac_threshold),
      min_triangulate_views_(min_triangulate_views),
      triangulate_error_threshold_(triangulate_error_threshold),
      number_of_threads_(number_of_threads) {}

  Err operator() (const ImageKeysetContainer& image_keysets,
                  const ImageIntrinsicMap& image_intrinsic_map,
                  const TrackContainer& tracks,
                  IntrinsicParamsContainer& intrinsic_params_set,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  ImageExtrinsicMap& image_extrinsic_map,
                  PointContainer& points,
                  TrackPointMap& track_point_map,
                  ViewInfoIndexer& view_info_indexer,
                  hs::progress::ProgressManager* progress_manager = NULL
                  ) const
  {
    return Run(image_keysets,
               image_intrinsic_map,
               tracks,
               intrinsic_params_set,
               extrinsic_params_set,
               image_extrinsic_map,
               points,
               track_point_map,
               view_info_indexer,
               progress_manager
               );
  }

  Err Run(const ImageKeysetContainer& image_keysets,
          const ImageIntrinsicMap& image_intrinsic_map,
          const TrackContainer& tracks,
          IntrinsicParamsContainer& intrinsic_params_set,
          ExtrinsicParamsContainer& extrinsic_params_set,
          ImageExtrinsicMap& image_extrinsic_map,
          PointContainer& points,
          TrackPointMap& track_point_map,
          ViewInfoIndexer& view_info_indexer,
          hs::progress::ProgressManager* progress_manager = NULL
          ) const
  {
    size_t number_of_images = image_keysets.size();
    if (image_intrinsic_map.Size() != number_of_images) return -1;
    ImageViewTracksContainer image_view_tracks_set;
    if (Initialize(number_of_images,
                   tracks,
                   image_view_tracks_set,
                   view_info_indexer) != 0)
    {
      return -1;
    }

    while (1)
    {
      if (progress_manager)	//若指定了ProgressManager
      {
        if (!progress_manager->CheckKeepWorking())
        {
          break;
        }
      }

      BundleAdjustmentOptimizor bundle_adjustment_optimizor(
                                  number_of_threads_);	//光束法平差校准
      if (bundle_adjustment_optimizor(image_keysets,
                                      image_intrinsic_map,
                                      tracks,
                                      image_extrinsic_map,
                                      track_point_map,
                                      view_info_indexer,
                                      intrinsic_params_set,
                                      extrinsic_params_set,
                                      points) != 0)
      {
        break;
      }

      ExtrinsicParamsContainer new_extrinsic_params_set;
      std::vector<size_t> new_image_ids;
      if (image_expandor_(image_keysets,
                          intrinsic_params_set,
                          image_intrinsic_map,
                          points,
                          track_point_map,
                          image_view_tracks_set,
                          image_extrinsic_map,
                          view_info_indexer,
                          new_extrinsic_params_set,
                          new_image_ids) != 0)	//调用ImageExpandor处理新照片的加入, 若无新照片则跳出
      {
        break;
      }

      for (size_t i = 0; i < new_extrinsic_params_set.size(); i++)
      {
        extrinsic_params_set.push_back(new_extrinsic_params_set[i]);
        image_extrinsic_map[new_image_ids[i]] = extrinsic_params_set.size() - 1;
      }

      //if (bundle_adjustment_optimizor(image_keysets,
      //                                image_intrinsic_map,
      //                                tracks,
      //                                image_extrinsic_map,
      //                                track_point_map,
      //                                view_info_indexer,
      //                                intrinsic_params_set,
      //                                extrinsic_params_set,
      //                                points) != 0)
      //{
      //  break;
      //}

      PointExpandor point_expandor;	//处理空间三维点的加入, 类似ImageExpandor
      if (point_expandor(image_keysets,
                         intrinsic_params_set,
                         image_intrinsic_map,
                         tracks,
                         extrinsic_params_set,
                         image_extrinsic_map,
                         min_triangulate_views_,
                         triangulate_error_threshold_,
                         points,
                         track_point_map,
                         view_info_indexer) != 0)
      {
        break;
      }

      if (progress_manager)
      {
        progress_manager->SetCurrentSubProgressCompleteRatio(
          float(extrinsic_params_set.size() / float(image_keysets.size())));
      }
    }

    return 0;
  }

protected:
  Err Initialize(size_t number_of_images,
                 const TrackContainer& tracks,
                 ImageViewTracksContainer& image_view_tracks_set,	//这个container是ImageViewTracks的向量, 而ImageViewTracks是vector<size_t>
                 ViewInfoIndexer& view_info_indexer) const
  {
    //计算每张影像拍到的track
    image_view_tracks_set.resize(number_of_images);
    size_t number_of_tracks = tracks.size();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();	//获取每个track的size
      for (size_t j = 0; j < number_of_views; j++)
      {
        image_view_tracks_set[tracks[i][j].first].push_back(i);
      }
    }

    return 0;
  }

protected:
  ImageExpandor image_expandor_;
  size_t min_triangulate_views_;
  Scalar triangulate_error_threshold_;
  size_t number_of_threads_;
};

}
}
}

#endif
