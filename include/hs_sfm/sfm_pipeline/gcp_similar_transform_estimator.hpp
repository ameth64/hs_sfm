#ifndef _HS_SFM_SFM_PIPELINE_GCP_SIMILAR_TRANSFORM_ESTIMATOR_HPP_
#define _HS_SFM_SFM_PIPELINE_GCP_SIMILAR_TRANSFORM_ESTIMATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/similar_transform_estimator.hpp"
#include "hs_sfm/sfm_pipeline/point_expandor.hpp"

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
class GCPSimilarTransformEstimator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

private:
  typedef hs::sfm::SimilarTransformEstimator<Scalar> SimilarTransformEstimator;

public:
  typedef hs::sfm::triangulate::MultipleViewMaximumLikelihoodEstimator<Scalar>
          TriangulateMLEstimator;
  typedef typename TriangulateMLEstimator::Key Key;
  typedef typename TriangulateMLEstimator::KeyContainer KeyContainer;
  typedef typename SimilarTransformEstimator::Rotation Rotation;
  typedef typename SimilarTransformEstimator::Translate Translate;

public:
  Err operator() (const KeysetContainer& keysets_gcp,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const ExtrinsicParamsContainer& extrinsic_params_set,
                  const hs::sfm::TrackContainer& tracks_gcp,
                  const PointContainer& gcps,
                  const hs::sfm::ObjectIndexMap& image_intrinsic_map,
                  const hs::sfm::ObjectIndexMap& image_extrinsic_map,
                  size_t min_triangulate_views,
                  Scalar triangulate_error_threshold,
                  Rotation& rotation_similar,
                  Translate& translate_similar,
                  Scalar& scale_similar,
                  hs::sfm::ObjectIndexMap& track_point_map_gcp,
                  PointContainer& gcps_relative) const
  {
    //计算像控点在相对坐标系下的坐标
    size_t number_of_gcps = tracks_gcp.size();
    if (gcps.size() != number_of_gcps)
    {
      std::cout<<"Error:-1\n";
      return -1;
    }
    track_point_map_gcp.Resize(number_of_gcps);
    gcps_relative.clear();
    for (size_t track_id = 0; track_id < tracks_gcp.size(); track_id++)
    {
      size_t number_of_views = tracks_gcp[track_id].size();
      hs::sfm::Track track_views;
      for (size_t i = 0; i < number_of_views; i++)
      {
        size_t image_id = tracks_gcp[track_id][i].first;
        if (image_extrinsic_map.IsValid(image_id))
        {
          track_views.push_back(tracks_gcp[track_id][i]);
        }
      }

      if (track_views.size() >= min_triangulate_views)
      {
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
          keys[i] = keysets_gcp[image_id][key_id];
        }

        TriangulateMLEstimator estimator;
        Point point;
        if (estimator(intrinsic_params_set_view,
                      extrinsic_params_set_view,
                      keys,
                      point) == 0)
        {
          gcps_relative.push_back(point);
          track_point_map_gcp[track_id] = gcps_relative.size() - 1;
        }
      }
    }

    //计算相似变换
    size_t number_of_available_gcps = gcps_relative.size();
    if (number_of_available_gcps < 4)
    {
      std::cout<<"Error:-3\n";
      return -1;
    }
    PointContainer gcps_absolute(number_of_available_gcps);
    for (size_t i = 0; i < number_of_gcps; i++)
    {
      if (track_point_map_gcp.IsValid(i))
      {
        gcps_absolute[track_point_map_gcp[i]] = gcps[i];
      }
    }

    SimilarTransformEstimator transform_estimator;
    if (transform_estimator(gcps_relative, gcps_absolute,
                            rotation_similar,
                            translate_similar,
                            scale_similar) != 0)
    {
      std::cout<<"Error:-4\n";
      return -1;
    }

#ifdef _DEBUG
    for (size_t i = 0; i < number_of_available_gcps; i++)
    {
      Point gcp_abs = gcps_relative[i];
      gcp_abs = scale_similar * (rotation_similar * gcp_abs) + translate_similar;
      Point diff = gcp_abs - gcps_absolute[i];

      int bp = 0;
    }
#endif

    return 0;
  }
};

}
}
}

#endif
