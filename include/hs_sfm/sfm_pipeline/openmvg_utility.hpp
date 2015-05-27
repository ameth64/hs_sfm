#ifndef _HS_SFM_SFM_PIPELINE_OPENMVG_UTILITY_HPP_
#define _HS_SFM_SFM_PIPELINE_OPENMVG_UTILITY_HPP_

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
class OpenMVGUtility
{
public:
  typedef _Scalar Scalar;

  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::math::geometry::Rotation3D<Scalar> Rotation;
  typedef EIGEN_STD_MAP(size_t, Rotation) RotationContainer;

public:
  int FillSFMData(const hs::sfm::ObjectIndexMap& image_intrinsic_map,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  openMVG::SfM_Data& sfm_data) const
  {
    for (size_t image_id = 0; image_id < image_intrinsic_map.Size();
         image_id++)
    {
      size_t intrinsic_id = image_intrinsic_map[image_id];
      const IntrinsicParams& intrinsic_params =
        intrinsic_params_set[intrinsic_id];
      std::shared_ptr<openMVG::View> view(new openMVG::View);
      view->id_view = openMVG::IndexT(image_id);
      view->id_intrinsic = openMVG::IndexT(intrinsic_id);
      view->id_pose = openMVG::IndexT(image_id);
      sfm_data.views.insert(std::make_pair(view->id_view, view));
    }
    for (size_t intrinsic_id = 0; intrinsic_id < intrinsic_params_set.size();
         intrinsic_id++)
    {
      const IntrinsicParams& intrinsic_params =
        intrinsic_params_set[intrinsic_id];
      std::shared_ptr<openMVG::Pinhole_Intrinsic> intrinsic(
        new openMVG::Pinhole_Intrinsic_Radial_K3(
          0, 0,
          double(intrinsic_params.focal_length()),
          double(intrinsic_params.principal_point_x()),
          double(intrinsic_params.principal_point_y()),
          double(intrinsic_params.k1()),
          double(intrinsic_params.k2()),
          double(intrinsic_params.k3())));
      sfm_data.intrinsics.insert(std::make_pair(openMVG::IndexT(intrinsic_id),
                                                intrinsic));
    }
    return 0;
  }

  int FillNormalizedFeaturesProvider(
    const KeysetContainer& keysets,
    const hs::sfm::ObjectIndexMap& image_intrinsic_map,
    const IntrinsicParamsContainer& intrinsic_params_set,
    openMVG::Features_Provider& normalized_features_provider) const
  {
    typedef typename IntrinsicParams::KMatrix KMatrix;
    typedef EIGEN_VECTOR(Scalar, 2) Key;
    typedef EIGEN_VECTOR(Scalar, 3) HKey;

    for (size_t image_id = 0; image_id < keysets.size(); image_id++)
    {
      openMVG::PointFeatures point_features;
      size_t intrinsic_id = image_intrinsic_map[image_id];
      const IntrinsicParams& intrinsic_params =
        intrinsic_params_set[intrinsic_id];
      KMatrix K_inverse = intrinsic_params.GetKMatrix().inverse();
      for (size_t key_id = 0; key_id < keysets[image_id].size(); key_id++)
      {
        Key key = keysets[image_id][key_id];
        HKey hkey;
        hkey << key(0),
                key(1),
                Scalar(1);
        hkey = K_inverse * hkey;
        hkey /= hkey[2];

        openMVG::PointFeature point_feature((float)(hkey[0]),
                                            (float)(hkey[1]));
        point_features.push_back(point_feature);
      }
      normalized_features_provider.feats_per_view.insert(
        std::make_pair(image_id, point_features));
    }
    return 0;
  }

  int FillFeaturesProvider(
    const KeysetContainer& keysets,
    openMVG::Features_Provider& features_provider) const
  {
    typedef EIGEN_VECTOR(Scalar, 2) Key;

    for (size_t image_id = 0; image_id < keysets.size(); image_id++)
    {
      openMVG::PointFeatures point_features;
      for (size_t key_id = 0; key_id < keysets[image_id].size(); key_id++)
      {
        Key key = keysets[image_id][key_id];
        openMVG::PointFeature point_feature((float)(key[0]),
                                            (float)(key[1]));
        point_features.push_back(point_feature);
      }
      features_provider.feats_per_view.insert(
        std::make_pair(image_id, point_features));
    }
    return 0;
  }

  int FillMatchesProvider(const MatchContainer& matches,
                          openMVG::Matches_Provider& matches_provider) const
  {
    auto itr_image_pair = matches.begin();
    auto itr_image_pair_end = matches.end();
    for (; itr_image_pair != itr_image_pair_end; ++itr_image_pair)
    {
      openMVG::Pair pair;
      pair.first = openMVG::IndexT(itr_image_pair->first.first);
      pair.second = openMVG::IndexT(itr_image_pair->first.second);

      openMVG::matching::IndMatches ind_matches;
      auto itr_key_pair = itr_image_pair->second.begin();
      auto itr_key_pair_end = itr_image_pair->second.end();
      for (; itr_key_pair != itr_key_pair_end; ++itr_key_pair)
      {
        openMVG::matching::IndMatch ind_match;
        ind_match._i = openMVG::IndexT(itr_key_pair->first);
        ind_match._j = openMVG::IndexT(itr_key_pair->second);
        ind_matches.push_back(ind_match);
      }
      matches_provider._pairWise_matches.insert(
        std::make_pair(pair, ind_matches));
    }
    return 0;
  }

  int FillMapGlobalR(
    const RotationContainer& global_rotations,
    openMVG::Hash_Map<openMVG::IndexT, openMVG::Mat3>& map_global_r) const
  {
    typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;
    auto itr_rotation = global_rotations.begin();
    auto itr_rotation_end = global_rotations.end();
    for (; itr_rotation != itr_rotation_end; ++itr_rotation)
    {
      RMatrix R = itr_rotation->second;
      map_global_r.insert(std::make_pair(openMVG::IndexT(itr_rotation->first),
                                         openMVG::Mat3(R)));
    }
    return 0;
  }
};

}
}
}

#endif
