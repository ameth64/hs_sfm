#ifndef _HS_SFM_SFM_PIPELINE_GLOBAL_SFM_HPP_
#define _HS_SFM_SFM_PIPELINE_GLOBAL_SFM_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_progress/progress_utility/progress_manager.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"
#include "hs_sfm/sfm_pipeline/best_pair_selector.hpp"
#include "hs_sfm/sfm_pipeline/scene_expandor.hpp"
#include "hs_sfm/sfm_pipeline/global_types.hpp"
#include "hs_sfm/sfm_pipeline/global_rotation_average.hpp"
#include "hs_sfm/sfm_pipeline/global_position_average.hpp"
#include "hs_sfm/sfm_pipeline/point_expandor.hpp"

#define TRY_1DSFM 1
#if TRY_1DSFM
#include <fstream>
#include <iostream>
#include <iomanip>
#include "hs_sfm/sfm_file_io/scene_ply_saver.hpp"
#endif

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
class GlobalSFM
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

private:
  typedef hs::sfm::pipeline::RotationPair<Scalar> RotationPair;
  typedef hs::sfm::pipeline::PositionPair<Scalar> PositionPair;
  typedef typename RotationPair::Index Index;
  struct EpipolarEdge
  {
    Index first_id;
    Index second_id;
    ExtrinsicParams extrinsic_params_relative;
  };
  typedef EIGEN_STD_VECTOR(EpipolarEdge) EpipolarEdgeContainer;

public:
  GlobalSFM()
    : number_of_threads_(1) {}

  GlobalSFM(size_t number_of_threads)
    : number_of_threads_(number_of_threads) {}

  Err operator() (const ObjectIndexMap& image_intrinsic_map,
                  const hs::sfm::MatchContainer& matches,
                  const KeysetContainer& keysets,
                  IntrinsicParamsContainer& intrinsic_params_set,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  hs::sfm::ObjectIndexMap& image_extrinsic_map,
                  PointContainer& points,
                  TrackContainer& tracks,
                  hs::sfm::ObjectIndexMap& track_point_map,
                  hs::sfm::ViewInfoIndexer& view_info_indexer,
                  hs::progress::ProgressManager* progress_manager = NULL) const
  {
    hs::sfm::MatchesTracksConvertor matches_track_convertor;
    if (matches_track_convertor(matches, tracks) != 0)
    {
      return -1;
    }
    auto itr_track = tracks.begin();
    auto itr_track_end = tracks.end();
    for (; itr_track != itr_track_end; ++itr_track)
    {
      std::sort(itr_track->begin(), itr_track->end());
    }
    std::sort(tracks.begin(), tracks.end());
    MatchContainer matches_filtered;
    if (matches_track_convertor(tracks, matches_filtered) != 0)
    {
      return -1;
    }

    size_t number_of_tracks = tracks.size();
    track_point_map.Resize(number_of_tracks);

    //构造view info indexer
    view_info_indexer.SetViewInfoByTracks(tracks);

    Err result = 0;

    while (1)
    {
      EpipolarEdgeContainer epipolar_edges;
      result = ComputeEpipolarEdges(image_intrinsic_map,
                                    matches,
                                    keysets,
                                    intrinsic_params_set,
                                    epipolar_edges);
      if (result != 0) break;

#if TRY_1DSFM


      std::cout<<"tracks.size():"<<tracks.size()<<"\n";
      //Output cc file.
      std::ofstream cc_file("cc.txt");
      for (size_t i = 0; i < image_intrinsic_map.Size(); i++)
      {
        cc_file<<i<<"\n";
      }
      cc_file.close();

      //Output eg file.
      std::ofstream eg_file("EGs.txt");
      auto itr_edge = epipolar_edges.begin();
      auto itr_edge_end = epipolar_edges.end();
      eg_file.setf(std::ios::fixed);
      eg_file<<std::setprecision(8);
      for (; itr_edge != itr_edge_end; ++itr_edge)
      {
        EIGEN_MATRIX(Scalar, 3, 3) R =
          itr_edge->extrinsic_params_relative.rotation();
        const Point& c = itr_edge->extrinsic_params_relative.position();
        Point t = -R * c;
        R.transposeInPlace();
        eg_file<<itr_edge->first_id<<" "<<itr_edge->second_id<<" "
               <<R(0, 0)<<" "<<R(0, 1)<<" "<<R(0, 2)<<" "
               <<R(1, 0)<<" "<<R(1, 1)<<" "<<R(1, 2)<<" "
               <<R(2, 0)<<" "<<R(2, 1)<<" "<<R(2, 2)<<" "
               <<c(0)<<" "<<c(1)<<" "<<c(2)<<"\n";
      }
      eg_file.close();

      //Output coords file.
      std::ofstream coords_file("coords.txt");
      coords_file.setf(std::ios::fixed);
      coords_file<<std::setprecision(6);
      for (size_t i = 0; i < image_intrinsic_map.Size(); i++)
      {
        size_t intrinsic_id = image_intrinsic_map[i];
        const IntrinsicParams& intrinsic_params =
          intrinsic_params_set[intrinsic_id];
        const Keyset& keyset = keysets[i];
        coords_file<<"#index = "<<i<<", name = "<<i<<".jpg, keys = "
                   <<keyset.size()<<", px = "
                   <<intrinsic_params.principal_point_x()
                   <<", py = "<<intrinsic_params.principal_point_y()
                   <<", focal = "<<intrinsic_params.focal_length()<<"\n";
        for (size_t j = 0; j < keyset.size(); j++)
        {
          coords_file<<j<<" "<<keyset[j][0]<<" "<<keyset[j][1]<<" 0 0 0 0 0\n";
        }

      }
      coords_file.close();

      //Output tracks file.
      std::ofstream tracks_file("tracks.txt");
      tracks_file<<tracks.size()<<"\n";
      for (size_t i = 0; i < tracks.size(); i++)
      {
        if (tracks[i].size() > 1)
        {
          tracks_file<<tracks[i].size();
          for (size_t j = 0; j < tracks[i].size(); j++)
          {
            tracks_file<<" "<<tracks[i][j].first<<" "<<tracks[i][j].second;
          }
          tracks_file<<"\n";
        }
      }
      tracks_file.close();

#define GETLINE(file) \
        std::getline(file, line);\
        ss.clear();\
        ss.str(line);

      std::stringstream ss;
      std::string line;

      //Input rot_solution file
      typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;
      EIGEN_STD_MAP(size_t, RMatrix) rotations;
      std::ifstream rot_file("rot_solution.txt");
      while (!rot_file.eof())
      {
        GETLINE(rot_file);
        if (line.empty()) break;

        size_t image_id;
        RMatrix R;
        ss>>image_id>>R(0, 0)>>R(0, 1)>>R(0, 2)
                    >>R(1, 0)>>R(1, 1)>>R(1, 2)
                    >>R(2, 0)>>R(2, 1)>>R(2, 2);
        R.transposeInPlace();
        R.row(1) *= Scalar(-1);
        R.row(2) *= Scalar(-1);
        if (rotations.find(image_id) == rotations.end())
        {
          rotations[image_id] = R;
        }
      }
      rot_file.close();

      //Input trans file
      EIGEN_STD_MAP(size_t, Point) positions;
      std::ifstream trans_file("trans_problem_solution.txt");
      while (!trans_file.eof())
      {
        GETLINE(trans_file);
        if (line.empty()) break;

        size_t image_id;
        Point position;
        ss>>image_id>>position(0)>>position(1)>>position(2);
        if ((image_id >= 0 && image_id < image_intrinsic_map.Size()) &&
            positions.find(image_id) == positions.end())
        {
          positions[image_id] = position;
        }
      }
      trans_file.close();

      image_extrinsic_map.Resize(image_intrinsic_map.Size());
      extrinsic_params_set.clear();
      for (size_t i = 0; i < image_extrinsic_map.Size(); i++)
      {
        auto itr_rotation = rotations.find(i);
        auto itr_position = positions.find(i);
        if (itr_rotation != rotations.end() && itr_position != positions.end())
        {
          ExtrinsicParams extrinsic_params;
          extrinsic_params.rotation() = itr_rotation->second;
          extrinsic_params.position() = itr_position->second;
          extrinsic_params_set.push_back(extrinsic_params);
          image_extrinsic_map[i] = extrinsic_params_set.size() - 1;
        }
      }

      typedef hs::sfm::fileio::ScenePLYSaver<Scalar, size_t> SceneSaver;
      typedef typename SceneSaver::Image Image;
      typedef typename SceneSaver::ImageContainer ImageContainer;
      SceneSaver scene_saver(0.5);
      size_t number_of_images = image_intrinsic_map.Size();
      IntrinsicParamsContainer intrinsic_params_set_saver(number_of_images);
      ExtrinsicParamsContainer extrinsic_params_set_saver(number_of_images);
      Image image;
      image.m_width = 6000;
      image.m_height = 4000;
      ImageContainer images(number_of_images, image);
      for (size_t i = 0; i < number_of_images; i++)
      {
        if (image_intrinsic_map.IsValid(i))
        {
          size_t intrinsic_id = image_intrinsic_map[i];
          intrinsic_params_set_saver[i] = intrinsic_params_set[intrinsic_id];
        }
        if (image_extrinsic_map.IsValid(i))
        {
          size_t extrinsic_id = image_extrinsic_map[i];
          extrinsic_params_set_saver[i] = extrinsic_params_set[extrinsic_id];
        }
      }

      scene_saver(std::string("global_scene.ply"),
                  intrinsic_params_set_saver,
                  extrinsic_params_set_saver,
                  images,
                  points);

#endif

      result = TriangluatePoints(keysets,
                                 intrinsic_params_set,
                                 image_intrinsic_map,
                                 tracks,
                                 extrinsic_params_set,
                                 image_extrinsic_map,
                                 points,
                                 track_point_map,
                                 view_info_indexer);
      if (result != 0) break;

#if TRY_1DSFM
      std::cout<<"points.size():"<<points.size()<<"\n";
#endif

      result = BundlerAdjustmentOptimize(keysets,
                                         image_intrinsic_map,
                                         tracks,
                                         image_extrinsic_map,
                                         track_point_map,
                                         view_info_indexer,
                                         intrinsic_params_set,
                                         extrinsic_params_set,
                                         points);
      if (result != 0) break;

      break;
    }

    return 0;
  }

private:
  Err ComputeEpipolarEdges(const ObjectIndexMap& image_intrinsic_map,
                           const hs::sfm::MatchContainer& matches,
                           const KeysetContainer& keysets,
                           const IntrinsicParamsContainer& intrinsic_params_set,
                           EpipolarEdgeContainer& epipolar_edges) const
  {
    typedef hs::sfm::essential::EMatrix5PointsRansacRefiner<Scalar>
            EMatrixRansacRefiner;
    typedef hs::sfm::essential::EMatrix5PointsCalculator<Scalar>
            EMatrixCalculator;
    typedef hs::sfm::essential::EMatrixExtrinsicParamsPointsCalculator<Scalar>
            ExtrinsicParamsPointsCalculator;
    typedef typename ExtrinsicParamsPointsCalculator::HKeyPair HKeyPair;
    typedef typename ExtrinsicParamsPointsCalculator::HKeyPairContainer
                     HKeyPairContainer;
    typedef typename ExtrinsicParamsPointsCalculator::EMatrix EMatrix;
    typedef typename IntrinsicParams::KMatrix KMatrix;
    typedef typename EMatrixRansacRefiner::IndexSet IndexSet;

    EMatrixRansacRefiner ransac_refiner;
    ExtrinsicParamsPointsCalculator extrinsic_params_calculator;
    auto itr_image_pair = matches.begin();
    auto itr_image_pair_end = matches.end();
    for (; itr_image_pair != itr_image_pair_end; ++itr_image_pair)
    {
      HKeyPairContainer key_pairs;
      size_t image_id_left = itr_image_pair->first.first;
      size_t image_id_right = itr_image_pair->first.second;
      size_t intrinsic_id_left = image_intrinsic_map[image_id_left];
      size_t intrinsic_id_right = image_intrinsic_map[image_id_right];
      const IntrinsicParams& intrinsic_params_left =
        intrinsic_params_set[intrinsic_id_left];
      const IntrinsicParams& intrinsic_params_right =
        intrinsic_params_set[intrinsic_id_right];
      KMatrix K_left_inverse = intrinsic_params_left.GetKMatrix().inverse();
      KMatrix K_right_inverse = intrinsic_params_right.GetKMatrix().inverse();
      auto itr_key_pair = itr_image_pair->second.begin();
      auto itr_key_pair_end = itr_image_pair->second.end();
      for (; itr_key_pair != itr_key_pair_end; ++itr_key_pair)
      {
        size_t key_left_id = itr_key_pair->first;
        size_t key_right_id = itr_key_pair->second;
        HKeyPair key_pair;
        key_pair.first.segment(0, 2) = keysets[image_id_left][key_left_id];
        key_pair.first[2] = Scalar(1);
        key_pair.second.segment(0, 2) = keysets[image_id_right][key_right_id];
        key_pair.second[2] = Scalar(1);
        key_pair.first = K_left_inverse * key_pair.first;
        key_pair.second = K_right_inverse * key_pair.second;
        key_pairs.push_back(key_pair);
      }

      EpipolarEdge epipolar_edge;
      epipolar_edge.first_id = image_id_left;
      epipolar_edge.second_id = image_id_right;
      HKeyPairContainer key_pairs_refined;
      IndexSet inlier_indices;
      EMatrix e_matrix;
      if (ransac_refiner(key_pairs, 8/ intrinsic_params_left.focal_length(),
                         key_pairs_refined, inlier_indices,
                         e_matrix) != 0) continue;

      PointContainer points_essential;
      if (extrinsic_params_calculator(
            e_matrix, key_pairs_refined,
            epipolar_edge.extrinsic_params_relative,
            points_essential) != 0) continue;

      epipolar_edges.push_back(epipolar_edge);
    }

    return 0;
  }

  Err RotationAverage() const
  {
    return 0;
  }

  Err TranslateAverage() const
  {
    return 0;
  }

  Err TriangluatePoints(
    const KeysetContainer& keysets,
    const IntrinsicParamsContainer& intrinsic_params_set,
    const hs::sfm::ObjectIndexMap& image_intrinsic_map,
    const hs::sfm::TrackContainer& tracks,
    const ExtrinsicParamsContainer& extrinsic_params_set,
    const hs::sfm::ObjectIndexMap& image_extrinsic_map,
    PointContainer& points,
    hs::sfm::ObjectIndexMap& track_point_map,
    hs::sfm::ViewInfoIndexer& view_info_indexer) const
  {
    typedef PointExpandor<Scalar> Expandor;
    Expandor expandor;
    return expandor(keysets,
                    intrinsic_params_set,
                    image_intrinsic_map,
                    tracks,
                    extrinsic_params_set,
                    image_extrinsic_map,
                    2, Scalar(64),
                    points,
                    track_point_map,
                    view_info_indexer);
  }

  Err BundlerAdjustmentOptimize(
    const KeysetContainer keysets,
    const hs::sfm::ObjectIndexMap& image_intrinsic_map,
    const TrackContainer& tracks,
    const hs::sfm::ObjectIndexMap& image_extrinsic_map,
    const hs::sfm::ObjectIndexMap& track_point_map,
    const ViewInfoIndexer& view_info_indexer,
    IntrinsicParamsContainer& intrinsic_params_set,
    ExtrinsicParamsContainer& extrinsic_params_set,
    PointContainer& points) const
  {
    typedef BundleAdjustmentOptimizor<Scalar> Optimizor;
    Optimizor optimizor(number_of_threads_);
    return optimizor(keysets,
                     image_intrinsic_map,
                     tracks,
                     image_extrinsic_map,
                     track_point_map,
                     view_info_indexer,
                     intrinsic_params_set,
                     extrinsic_params_set,
                     points);
  }

private:
  size_t number_of_threads_;
};

}
}
}

#endif
