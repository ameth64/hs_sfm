#ifndef _HS_SFM_INCREMENTAL_POINT_EXPANDOR_HPP_
#define _HS_SFM_INCREMENTAL_POINT_EXPANDOR_HPP_

#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/projective_functions.hpp"
#include "hs_sfm/triangulate/multiple_view_maximum_likelihood_estimator.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar>
class PointExpandor
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

private:
  typedef hs::sfm::triangulate::MultipleViewMaximumLikelihoodEstimator<Scalar>
          TriangulateMLEstimator;
  typedef typename TriangulateMLEstimator::Key Key;
  typedef typename TriangulateMLEstimator::KeyContainer KeyContainer;
  typedef hs::sfm::ProjectiveFunctions<Scalar> ProjectiveFunctionsType;

public:
  Err operator() (const ImageKeysetContainer& image_keysets,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const ImageIntrinsicMap& image_intrinsic_map,
                  const TrackContainer& tracks,
                  const ExtrinsicParamsContainer& extrinsic_params_set,
                  const ImageExtrinsicMap& image_extrinsic_map,
                  size_t min_triangulate_views,
                  Scalar triangulate_error_threshold,
                  PointContainer& points,
                  TrackPointMap& track_point_map,
                  ViewInfoIndexer& view_info_indexer) const
  {
    //获取仍未被加入的三维点
    size_t number_of_tracks = tracks.size();

    for (size_t track_id = 0; track_id < number_of_tracks; track_id++)
    {
      if (!track_point_map.IsValid(track_id))
      {
        size_t number_of_views = tracks[track_id].size();

        hs::sfm::Track track_views;
        for (size_t i = 0; i < number_of_views; i++)
        {
          size_t image_id = tracks[track_id][i].first;
          const ViewInfo& view_info =
            view_info_indexer.GetViewInfoByTrackImage(track_id, image_id);
          if (image_extrinsic_map.IsValid(image_id) &&
              !view_info.is_blunder)
          {
            track_views.push_back(tracks[track_id][i]);
          }
        }

        if (track_views.size() >= min_triangulate_views)
        {
          //三角化计算三维点
          size_t track_size = track_views.size();
          IntrinsicParamsContainer intrinsic_params_set_view(track_size);
          ExtrinsicParamsContainer extrinsic_params_set_view(track_size);
          KeyContainer keys(track_size);
          for (size_t i = 0; i < track_size; i++)
          {
            size_t image_id = track_views[i].first;
            size_t key_id = track_views[i].second;
            size_t intrinsic_id = image_intrinsic_map[image_id];
            size_t extrinsic_id = image_extrinsic_map[image_id];
            intrinsic_params_set_view[i] =
              intrinsic_params_set[intrinsic_id];
            extrinsic_params_set_view[i] =
              extrinsic_params_set[extrinsic_id];
            keys[i] = image_keysets[image_id][key_id];
          }

          TriangulateMLEstimator estimator;
          Point point;
          if (estimator(intrinsic_params_set_view,
                        extrinsic_params_set_view,
                        keys,
                        point) == 0)
          {
            //计算重投影误差
            bool is_blunder = false;
            for (size_t i = 0; i < track_size; i++)
            {
              Key key = ProjectiveFunctionsType::WorldPointProjectToImageKey(
                intrinsic_params_set_view[i], extrinsic_params_set_view[i],
                point);

              Scalar error = (keys[i] - key).norm();

              if (error > triangulate_error_threshold)
              {
                is_blunder = true;
                break;
              }
            }//for (size_t i = 0; i < track_size; i++)

            if (!is_blunder)
            {
              points.push_back(point);
              track_point_map[track_id] = points.size() - 1;
            }//if (!is_blunder)
            else
            {
              for (size_t i = 0; i < number_of_views; i++)
              {
                size_t image_id = tracks[track_id][i].first;
                ViewInfo& view_info =
                  view_info_indexer.GetViewInfoByTrackImage(track_id,
                  image_id);
                view_info.is_blunder = true;
              }
            }
          }
        }//if (track_views.size() >= min_triangulate_views_)
      }//if (!track_point_map.IsValid(track_id))
    }//for (size_t track_id = 0; track_id < number_of_tracks; track_id++)

    return 0;
  }

};

}
}
}

#endif
