#ifndef _HS_SFM_INCREMENTAL_GCP_SIMILAR_TRANSFORM_ESTIMATOR_HPP_
#define _HS_SFM_INCREMENTAL_GCP_SIMILAR_TRANSFORM_ESTIMATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/similar_transform_estimator.hpp"
#include "hs_sfm/incremental/point_expandor.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
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
  typedef hs::sfm::incremental::PointExpandor<Scalar> PointTriangulator;

public:
  typedef typename SimilarTransformEstimator::Rotation Rotation;
  typedef typename SimilarTransformEstimator::Translate Translate;

public:
  Err operator() (const KeysetContainer& keysets_gcp,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const ExtrinsicParamsContainer& extrinsic_params_set,
                  const hs::sfm::TrackContainer& tracks_gcp,
                  const PointContainer& gcps,
                  const hs::sfm::ObjectIndexMap& image_extrinsic_map,
                  size_t min_triangulate_views,
                  Scalar triangulate_error_threshold,
                  Rotation& rotation_similar,
                  Translate& translate_similar,
                  Scalar& scale_similar,
                  hs::sfm::ObjectIndexMap& track_point_map_gcp,
                  hs::sfm::ViewInfoIndexer& view_info_indexer_gcp,
                  PointContainer& gcps_relative) const
  {
    //计算像控点在相对坐标系下的坐标
    size_t number_of_gcps = tracks_gcp.size();
    if (gcps.size() != number_of_gcps)
    {
      return -1;
    }
    track_point_map_gcp.Resize(number_of_gcps);
    PointTriangulator triangulator;
    gcps_relative.clear();
    if (triangulator(keysets_gcp,
                     intrinsic_params_set,
                     tracks_gcp,
                     extrinsic_params_set,
                     image_extrinsic_map,
                     min_triangulate_views,
                     triangulate_error_threshold,
                     gcps_relative,
                     track_point_map_gcp,
                     view_info_indexer_gcp) != 0)
    {
      return -1;
    }

    //计算相似变换
    size_t number_of_available_gcps = gcps_relative.size();
    if (number_of_available_gcps < 4)
    {
      return -1;
    }
    PointContainer gcps_absolute;
    for (size_t i = 0; i < number_of_gcps; i++)
    {
      if (track_point_map_gcp.IsValid(i))
      {
        gcps_absolute.push_back(gcps[i]);
      }
    }

    SimilarTransformEstimator transform_estimator;
    if (transform_estimator(gcps_relative, gcps_absolute,
                            rotation_similar,
                            translate_similar,
                            scale_similar) != 0)
    {
      return -1;
    }

    return 0;
  }
};

}
}
}

#endif
