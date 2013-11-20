#ifndef _HS_SFM_INCREMENTAL_SCENE_EXPANDOR_HPP_
#define _HS_SFM_INCREMENTAL_SCENE_EXPANDOR_HPP_

#include <vector>
#include <map>

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
    Scalar triangulate_error_threshold = 4.0)
    : image_expandor_(add_new_image_matches_threshold,
                      pmatrix_ransac_threshold),
      min_triangulate_views_(min_triangulate_views),
      triangulate_error_threshold_(triangulate_error_threshold){}

  Err operator() (const ImageKeysetContainer& image_keysets,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const TrackContainer& tracks,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  ImageExtrinsicMap& image_extrinsic_map,
                  PointContainer& points,
                  TrackPointMap& track_point_map,
                  ViewInfoIndexer& view_info_indexer) const
  {
    return Run(image_keysets,
               intrinsic_params_set,
               tracks,
               extrinsic_params_set,
               image_extrinsic_map,
               points,
               track_point_map,
               view_info_indexer);
  }

  Err Run(const ImageKeysetContainer& image_keysets,
          const IntrinsicParamsContainer& intrinsic_params_set,
          const TrackContainer& tracks,
          ExtrinsicParamsContainer& extrinsic_params_set,
          ImageExtrinsicMap& image_extrinsic_map,
          PointContainer& points,
          TrackPointMap& track_point_map,
          ViewInfoIndexer& view_info_indexer) const
  {
    ImageViewTracksContainer image_view_tracks_set;
    if (Initialize(image_keysets,
                   intrinsic_params_set,
                   tracks,
                   image_extrinsic_map,
                   track_point_map,
                   image_view_tracks_set,
                   view_info_indexer) != 0)
    {
      return -1;
    }

    while (1)
    {
      BundleAdjustmentOptimizor bundle_adjustment_optimizor;
      if (bundle_adjustment_optimizor(image_keysets,
                                      intrinsic_params_set,
                                      tracks,
                                      image_extrinsic_map,
                                      track_point_map,
                                      view_info_indexer,
                                      extrinsic_params_set,
                                      points) != 0)
      {
        break;
      }

      ExtrinsicParams new_extrinsic_params;
      size_t new_image_id;
      if (image_expandor_(image_keysets,
                          intrinsic_params_set,
                          points,
                          track_point_map,
                          image_view_tracks_set,
                          image_extrinsic_map,
                          view_info_indexer,
                          new_extrinsic_params,
                          new_image_id) != 0)
      {
        break;
      }

      extrinsic_params_set.push_back(new_extrinsic_params);
      image_extrinsic_map[new_image_id] = extrinsic_params_set.size() - 1;

      PointExpandor point_expandor;
      if (point_expandor(image_keysets,
                         intrinsic_params_set,
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
    }

    return 0;
  }

protected:
  Err Initialize(const ImageKeysetContainer& image_keysets,
                 const IntrinsicParamsContainer& intrinsic_params_set,
                 const TrackContainer& tracks,
                 const ImageExtrinsicMap& image_extrinsic_map,
                 const TrackPointMap& track_point_map,
                 ImageViewTracksContainer& image_view_tracks_set,
                 ViewInfoIndexer& view_info_indexer) const
  {
    size_t number_of_images = image_keysets.size();
    if (number_of_images != intrinsic_params_set.size())
    {
      return -1;
    }

    //计算每个内参数对应的影像拍到的track
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

    //构造view info indexer
    view_info_indexer.SetViewInfoByTracks(tracks);

    return 0;
  }

protected:
  ImageExpandor image_expandor_;
  size_t min_triangulate_views_;
  Scalar triangulate_error_threshold_;
};

}
}
}

#endif
