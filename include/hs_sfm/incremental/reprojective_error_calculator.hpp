#ifndef _HS_SFM_INCREMENTAL_REPROJECTIVE_ERROR_CALCULATOR_HPP_
#define _HS_SFM_INCREMENTAL_REPROJECTIVE_ERROR_CALCULATOR_HPP_

#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar>
class ReprojectiveErrorCalculator
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

private:
  typedef CameraFunctions<Scalar> CameraFunctions;
  typedef typename CameraFunctions::ProjectionMatrix PMatrix;

public:
  Scalar operator() (const ImageKeysetContainer& image_keysets,
                     const IntrinsicParamsContainer& intrinsic_params_set,
                     const TrackContainer& tracks,
                     const ImageExtrinsicMap& image_extrinsic_map,
                     const TrackPointMap& track_point_map,
                     const ViewInfoIndexer& view_info_indexer,
                     const ExtrinsicParamsContainer& extrinsic_params_set,
                     const PointContainer& points) const
  {
    size_t number_of_tracks = tracks.size();
    size_t number_of_reprojections = 0;
    Scalar mean_reprojection_error = Scalar(0);
    for (size_t track_id = 0; track_id < number_of_tracks; track_id++)
    {
      if (track_point_map.IsValid(track_id))
      {
        size_t number_of_views = tracks[track_id].size();
        hs::sfm::Track track_bundle;
        for (size_t i = 0; i < number_of_views; i++)
        {
          size_t image_id = tracks[track_id][i].first;
          const ViewInfo& view_info =
            view_info_indexer.GetViewInfoByTrackImage(track_id, image_id);
          if (image_extrinsic_map.IsValid(image_id) &&
            !view_info.is_blunder)
          {
            size_t key_id = tracks[track_id][i].second;
            size_t exrinsic_id = image_extrinsic_map[image_id];
            size_t point_id = track_point_map[track_id];
            const ExtrinsicParams& extrinsic_params =
              extrinsic_params_set[exrinsic_id];
            const IntrinsicParams& intrinsic_params =
              intrinsic_params_set[image_id];
            EIGEN_VECTOR(Scalar, 2) key = image_keysets[image_id][key_id];
            PMatrix P = CameraFunctions::GetProjectionMatrix(intrinsic_params,
                                                             extrinsic_params);
            const Point& point = points[point_id];
            EIGEN_VECTOR(Scalar, 3) hkey = P.block(0, 0, 3, 3) * point +
                                           P.col(3);
            hkey /= hkey(2);
            Scalar reprojection_error = (key - hkey.segment(0, 2)).norm();
            mean_reprojection_error += reprojection_error;
            number_of_reprojections++;
          }
        }
      }
    }

    mean_reprojection_error /= Scalar(number_of_reprojections);

    return mean_reprojection_error;
  }
};

}
}
}

#endif
