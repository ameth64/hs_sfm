#ifndef _HS_SFM_SFM_PIPELINE_EPIPOLAR_EDGES_CALCULATOR_OPENMVG_HPP_
#define _HS_SFM_SFM_PIPELINE_EPIPOLAR_EDGES_CALCULATOR_OPENMVG_HPP_

#include "openMVG/types.hpp"
#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/essential.hpp"
#include "third_party/progress/progress.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_pipeline/global_types.hpp"
#include "hs_sfm/sfm_pipeline/openmvg_utility.hpp"

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
class EpipolarEdgesCalculatorOpenMVG
{
public:
  typedef _Scalar Scalar;
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::sfm::pipeline::EpipolarEdge<Scalar> EpipolarEdge;
  typedef EIGEN_STD_VECTOR(EpipolarEdge) EpipolarEdgeContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;
public:
  int operator() (const ObjectIndexMap& image_intrinsic_map,
                  const hs::sfm::MatchContainer& matches,
                  const KeysetContainer& keysets,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  EpipolarEdgeContainer& epipolar_edges) const
  {
    typedef hs::sfm::pipeline::OpenMVGUtility<Scalar> OpenMVGUtility;
    typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;

    OpenMVGUtility openmvg_utility;

    openMVG::SfM_Data sfm_data;
    openmvg_utility.FillSFMData(image_intrinsic_map, intrinsic_params_set,
                                sfm_data);

    openMVG::Features_Provider normalized_features_provider;
    openmvg_utility.FillNormalizedFeaturesProvider(
      keysets, image_intrinsic_map, intrinsic_params_set,
      normalized_features_provider);

    openMVG::Features_Provider features_provider;
    openmvg_utility.FillFeaturesProvider(
      keysets, features_provider);

    openMVG::Matches_Provider matches_provider;
    openmvg_utility.FillMatchesProvider(matches, matches_provider);

    //compute max connected component.
    const openMVG::Pair_Set pairs = matches_provider.getPairs();
    const std::set<openMVG::IndexT> set_remaining_ids =
      openMVG::graphUtils::CleanGraph_KeepLargestBiEdge_Nodes<
        openMVG::Pair_Set, openMVG::IndexT>(pairs);

    if (set_remaining_ids.empty())
    {
      return -1;
    }

    openMVG::KeepOnlyReferencedElement(set_remaining_ids,
                                       matches_provider._pairWise_matches);

    //compute relative pose.
    const openMVG::Pair_Set& pair_set = matches_provider.getPairs();
    const openMVG::Pair_Vec pair_vec(pair_set.begin(), pair_set.end());

    C_Progress_display my_progress_bar(
      pair_vec.size(), std::cout, "\nCompute_Relative_Rotations\n" );
    for (int i = 0; i < pair_vec.size(); ++i)
    {
      ++my_progress_bar;
      const openMVG::Pair current_pair = pair_vec[i];
      const size_t I = current_pair.first;
      const size_t J = current_pair.second;

      const openMVG::View * view_I = sfm_data.views[I].get();
      const openMVG::View * view_J = sfm_data.views[J].get();

      // Check that valid camera are existing for the pair index
      if (sfm_data.getIntrinsics().find(view_I->id_intrinsic) ==
          sfm_data.getIntrinsics().end() ||
          sfm_data.getIntrinsics().find(view_J->id_intrinsic) ==
          sfm_data.getIntrinsics().end())
        continue;

      const openMVG::IndMatches & vec_matchesInd =
        matches_provider._pairWise_matches.find(current_pair)->second;

      openMVG::Mat x1(2, vec_matchesInd.size()), x2(2, vec_matchesInd.size());
      for (size_t k = 0; k < vec_matchesInd.size(); ++k)
      {
        x1.col(k) =
          normalized_features_provider.
            feats_per_view[I][vec_matchesInd[k]._i].coords().cast<double>();
        x2.col(k) =
          normalized_features_provider.
            feats_per_view[J][vec_matchesInd[k]._j].coords().cast<double>();
      }

      const openMVG::IntrinsicBase * cam_I =
        sfm_data.getIntrinsics().at(view_I->id_intrinsic).get();
      const openMVG::IntrinsicBase * cam_J =
        sfm_data.getIntrinsics().at(view_J->id_intrinsic).get();
      if ( !isValid(cam_I->getType()) || !isValid(cam_J->getType()))
      {
        continue;
      }

      // Compute relative poses from the view graph
      // (thanks to a robust essential matrix estimation):

      openMVG::RelativePose_Info relativePose_info;
      // Compute max authorized error as geometric mean of camera plane
      // tolerated residual error
      relativePose_info.initial_residual_tolerance = std::pow(
        cam_I->imagePlane_toCameraPlaneError(2.5) +
        cam_J->imagePlane_toCameraPlaneError(2.5),
        1./2.);
      // Since we use normalized features:
      const std::pair<size_t, size_t> imageSize_I(1., 1.), imageSize_J(1.,1.);
      const openMVG::Mat3 K  = openMVG::Mat3::Identity();

      if (!robustRelativePose(K, K, x1, x2, relativePose_info,
                              imageSize_I, imageSize_J, 256))
      {
        continue;
      }
      bool bRefine_using_BA = true;
      if (bRefine_using_BA)
      {
        // Refine the defined scene
        openMVG::SfM_Data tiny_scene;
        tiny_scene.views.insert(*sfm_data.getViews().find(view_I->id_view));
        tiny_scene.views.insert(*sfm_data.getViews().find(view_J->id_view));
        tiny_scene.intrinsics.insert(
          *sfm_data.getIntrinsics().find(view_I->id_intrinsic));
        tiny_scene.intrinsics.insert(
          *sfm_data.getIntrinsics().find(view_J->id_intrinsic));

        // Init poses
        const openMVG::Pose3 & Pose_I =
          tiny_scene.poses[view_I->id_pose] =
            Pose3(openMVG::Mat3::Identity(), openMVG::Vec3::Zero());
        const openMVG::Pose3 & Pose_J =
          tiny_scene.poses[view_J->id_pose] =
            relativePose_info.relativePose;

        // Init structure
        const openMVG::Mat34 P1 = cam_I->get_projective_equivalent(Pose_I);
        const openMVG::Mat34 P2 = cam_J->get_projective_equivalent(Pose_J);
        openMVG::Landmarks & landmarks = tiny_scene.structure;
        for (size_t k = 0; k < x1.cols(); ++k) {
          const openMVG::Vec2 x1_ =
            features_provider.
              feats_per_view[I][vec_matchesInd[k]._i].coords().cast<double>();
          const openMVG::Vec2 x2_ =
            features_provider.
              feats_per_view[J][vec_matchesInd[k]._j].coords().cast<double>();
          openMVG::Vec3 X;
          openMVG::TriangulateDLT(P1, x1_, P2, x2_, &X);
          openMVG::Observations obs;
          obs[view_I->id_view] = Observation(x1_, vec_matchesInd[k]._i);
          obs[view_J->id_view] = Observation(x2_, vec_matchesInd[k]._j);
          landmarks[k].obs = obs;
          landmarks[k].X = X;
        }
        // - refine only Structure and
        // Rotations & translations (keep intrinsic constant)
        openMVG::Bundle_Adjustment_Ceres::BA_options options(false, false);
        options._linear_solver_type = ceres::DENSE_SCHUR;
        openMVG::Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        if (bundle_adjustment_obj.Adjust(tiny_scene, true, true, false))
        {
          // --> to debug: save relative pair geometry on disk
          // std::ostringstream os;
          // os << current_pair.first << "_" << current_pair.second << ".ply";
          // Save(tiny_scene, os.str(), ESfM_Data(STRUCTURE | EXTRINSICS));
          //
          const openMVG::Mat3 R1 =
            tiny_scene.poses[view_I->id_pose].rotation();
          const openMVG::Mat3 R2 =
            tiny_scene.poses[view_J->id_pose].rotation();
          const openMVG::Vec3 t1 =
            tiny_scene.poses[view_I->id_pose].translation();
          const openMVG::Vec3 t2 =
            tiny_scene.poses[view_J->id_pose].translation();
          // Compute relative motion and save it
          openMVG::Mat3 Rrel;
          openMVG::Vec3 trel;
          openMVG::RelativeCameraMotion(R1, t1, R2, t2, &Rrel, &trel);
          // Update found relative pose
          relativePose_info.relativePose =
            openMVG::Pose3(Rrel, -Rrel.transpose() * trel);
        }
      }

      {
        EpipolarEdge edge;
        edge.first_id = size_t(current_pair.first);
        edge.second_id = size_t(current_pair.second);
        RMatrix R = relativePose_info.relativePose.rotation();
        edge.extrinsic_params_relative.rotation() = R;
        edge.extrinsic_params_relative.position() =
          relativePose_info.relativePose.center();
        epipolar_edges.push_back(edge);
      }
    }

    return 0;
  }
};

}
}
}

#endif
