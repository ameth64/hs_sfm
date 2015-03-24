#ifndef _HS_SFM_INCREMENTAL_SCENE_EXPANDOR_HPP_
#define _HS_SFM_INCREMENTAL_SCENE_EXPANDOR_HPP_

#include <vector>
#include <map>

#include "hs_progress/progress_utility/progress_manager.hpp"

#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/incremental/image_expandor.hpp"
#include "hs_sfm/incremental/point_expandor.hpp"
#include "hs_sfm/incremental/bundle_adjustment_optimizor.hpp"

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

  typedef ObjectIndexMap TrackPointMap;
  typedef ObjectIndexMap ImageIntrinsicMap;
  typedef ObjectIndexMap ImageExtrinsicMap;

protected:
  typedef ImageExpandor<Scalar> ImageExpandor;
  typedef typename ImageExpandor::ImageViewTracks ImageViewTracks;
  typedef typename ImageExpandor::ImageViewTracksContainer
                   ImageViewTracksContainer;

  typedef PointExpandor<Scalar> PointExpandor;

  typedef BundleAdjustmentOptimizor<Scalar> BundleAdjustmentOptimizor;

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
                  hs::progress::ProgressManager* progress_manager = NULL) const
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
               progress_manager);
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
          hs::progress::ProgressManager* progress_manager = NULL) const
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
      if (progress_manager)
      {
        if (!progress_manager->CheckKeepWorking())
        {
          break;
        }
      }

      std::chrono::time_point<std::chrono::system_clock> start, end;
      start = std::chrono::system_clock::now();
      BundleAdjustmentOptimizor bundle_adjustment_optimizor(
                                  number_of_threads_);
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
      end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      std::cout<<"Bundle Adjustment took "<<elapsed_seconds.count()
               <<" seconds.\n";

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
                          new_image_ids) != 0)
      {
        break;
      }

      for (size_t i = 0; i < new_extrinsic_params_set.size(); i++)
      {
        std::cout<<"Added Image "<<new_image_ids[i]<<"\n";
        extrinsic_params_set.push_back(new_extrinsic_params_set[i]);
        image_extrinsic_map[new_image_ids[i]] = extrinsic_params_set.size() - 1;
      }
      std::cout<<new_image_ids.size()<<" images added.\n";

      PointExpandor point_expandor;
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
                 ImageViewTracksContainer& image_view_tracks_set,
                 ViewInfoIndexer& view_info_indexer) const
  {
    //计算每张影像拍到的track
    image_view_tracks_set.resize(number_of_images);
    size_t number_of_tracks = tracks.size();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();
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
