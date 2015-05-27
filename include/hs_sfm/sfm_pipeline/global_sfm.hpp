#ifndef _HS_SFM_SFM_PIPELINE_GLOBAL_SFM_HPP_
#define _HS_SFM_SFM_PIPELINE_GLOBAL_SFM_HPP_

#include <algorithm>

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
//#include "hs_sfm/sfm_pipeline/epipolar_edges_calculator.hpp"
#include "hs_sfm/sfm_pipeline/epipolar_edges_calculator_openmvg.hpp"

#define DEBUG_TMP 1
#if DEBUG_TMP
#include <iostream>
#endif

#define TRY_1DSFM 0
#if TRY_1DSFM
#include <fstream>
#include <iostream>
#include <iomanip>
#include "hs_sfm/sfm_file_io/scene_ply_saver.hpp"
#include "hs_sfm/sfm_utility/projective_functions.hpp"
#include "hs_sfm/fundamental/linear_8_points_ransac_refiner.hpp"
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
  typedef hs::sfm::pipeline::EpipolarEdge<Scalar> EpipolarEdge;
  typedef EpipolarEdgesCalculatorOpenMVG<Scalar> EdgesCalculator;
  typedef typename EdgesCalculator::EpipolarEdgeContainer EpipolarEdgeContainer;
  typedef typename RotationPair::Index Index;
#if TRY_1DSFM
  typedef typename Keyset::Key Key;
  typedef typename EIGEN_VECTOR(Scalar, 3) HKey;
  typedef hs::sfm::ProjectiveFunctions<Scalar> ProjectiveFunctions;
  typedef hs::sfm::fundamental::Linear8PointsRansacRefiner<Scalar> FRefiner;
  typedef typename FRefiner::KeyPairContainer FKeyPairContainer;
  typedef typename FRefiner::KeyPair FKeyPair;
#endif

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
        Point c = itr_edge->extrinsic_params_relative.position();
        R.transposeInPlace();
        R.col(1) *= Scalar(-1);
        R.col(2) *= Scalar(-1);
        R.row(1) *= Scalar(-1);
        R.row(2) *= Scalar(-1);
        c[1] *= Scalar(-1);
        c[2] *= Scalar(-1);
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
      SceneSaver scene_saver(0.1);
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

      image_extrinsic_map.Resize(image_intrinsic_map.Size());
      result = RotationAverage(matches,
                               epipolar_edges,
                               image_extrinsic_map,
                               extrinsic_params_set);
      if (result != 0) break;

      result = TranslateAverage(keysets,
                                image_intrinsic_map,
                                intrinsic_params_set,
                                matches,
                                epipolar_edges,
                                image_extrinsic_map,
                                extrinsic_params_set);
      if (result != 0) break;

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
    EdgesCalculator calculator;

    return calculator(image_intrinsic_map,
                      matches,
                      keysets,
                      intrinsic_params_set,
                      epipolar_edges);
  }

  Err RotationAverage(const hs::sfm::MatchContainer& matches,
                      const EpipolarEdgeContainer& epipolar_edges,
                      hs::sfm::ObjectIndexMap& image_extrinsic_map,
                      ExtrinsicParamsContainer& extrinsic_params_set) const
  {
    typedef GlobalRotationAverage<Scalar> RotationAverage;
    typedef typename RotationAverage::RotationPairContainer
            RotationPairContainer;
    typedef typename RotationAverage::IndexContainer IndexContainer;
    typedef typename RotationAverage::RotationContainer RotationContainer;

#if DEBUG_TMP
    std::cout<<"RotationAverage Start\n";
#endif

    RotationPairContainer rotation_pairs(epipolar_edges.size());
    std::vector<size_t> matches_counts;
    matches_counts.reserve(epipolar_edges.size());
    for (size_t i = 0; i < epipolar_edges.size(); i++)
    {
      rotation_pairs[i].first_id = epipolar_edges[i].first_id;
      rotation_pairs[i].second_id = epipolar_edges[i].second_id;
      rotation_pairs[i].rotation =
        epipolar_edges[i].extrinsic_params_relative.rotation();
      auto itr_image_pair =
        matches.find(std::make_pair(epipolar_edges[i].first_id,
                                    epipolar_edges[i].second_id));
      if (itr_image_pair != matches.end())
      {
        matches_counts.push_back(itr_image_pair->second.size());
      }
    }

#if DEBUG_TMP
    std::cout<<"epipolar_edges.size():"<<epipolar_edges.size()<<"\n";
    std::cout<<"Finish fill rotation pairs.\n";
#endif

    std::partial_sort(matches_counts.begin(),
                      matches_counts.begin() + matches_counts.size() / 2,
                      matches_counts.end());
    size_t median_count = matches_counts[matches_counts.size() / 2 - 1];

#if DEBUG_TMP
    std::cout<<"median_count:"<<median_count<<"\n";
#endif

    for (size_t i = 0; i < rotation_pairs.size(); i++)
    {
      auto itr_image_pair =
        matches.find(std::make_pair(epipolar_edges[i].first_id,
                                    epipolar_edges[i].second_id));
      if (itr_image_pair != matches.end())
      {
        rotation_pairs[i].weight =
          Scalar(itr_image_pair->second.size()) / Scalar(median_count);
      }
      else
      {
        rotation_pairs[i].weight = Scalar(1);
      }
    }

#if DEBUG_TMP
    std::cout<<"Finish compute weight.\n";
#endif

    RotationAverage rotation_average;
    IndexContainer global_rotation_indices;
    RotationContainer global_rotations;
    Err result = rotation_average(rotation_pairs,
                                  global_rotation_indices,
                                  global_rotations);
#if DEBUG_TMP
    std::cout<<"Finish rotation average.\n";
#endif
    if (result != 0) return result;

#if DEBUG_TMP
    std::cout<<"Copy rotation average output.\n";
#endif

    for (size_t i = 0; i < global_rotation_indices.size(); i++)
    {
      size_t image_id = global_rotation_indices[i];
      image_extrinsic_map[image_id] = i;
      ExtrinsicParams extrinsic_params;
      extrinsic_params.rotation() = global_rotations[i];
      extrinsic_params_set.push_back(extrinsic_params);
    }

#if DEBUG_TMP
    std::cout<<"RotationAverage End\n";
#endif
    return 0;
  }

  Err TranslateAverage(const KeysetContainer& keysets,
                       const hs::sfm::ObjectIndexMap& image_intrinsic_map,
                       const IntrinsicParamsContainer& intrinsic_params_set,
                       const hs::sfm::MatchContainer& matches,
                       const EpipolarEdgeContainer& epipolar_edges,
                       hs::sfm::ObjectIndexMap& image_extrinsic_map,
                       ExtrinsicParamsContainer& extrinsic_params_set) const
  {
    typedef GlobalPositionAverage<Scalar> PositionAverage;
    typedef typename PositionAverage::PositionPairContainer
            PositionPairContainer;
    typedef typename PositionAverage::IndexContainer IndexContainer;
    typedef typename PositionAverage::PositionContainer PositionContainer;
    typedef typename PositionAverage::RotationContainer RotationContainer;

#if DEBUG_TMP
    std::cout<<"TranslateAverage Start\n";
#endif

    PositionPairContainer position_pairs;
    std::vector<size_t> matches_counts;
    for (size_t i = 0; i < epipolar_edges.size(); i++)
    {
      if (image_extrinsic_map.IsValid(epipolar_edges[i].first_id) &&
          image_extrinsic_map.IsValid(epipolar_edges[i].second_id))
      {
        PositionPair position_pair;
        position_pair.first_id = epipolar_edges[i].first_id;
        position_pair.second_id = epipolar_edges[i].second_id;
        position_pair.position =
          epipolar_edges[i].extrinsic_params_relative.position();
        auto itr_image_pair =
          matches.find(std::make_pair(epipolar_edges[i].first_id,
                                      epipolar_edges[i].second_id));
        if (itr_image_pair != matches.end())
        {
          matches_counts.push_back(itr_image_pair->second.size());
        }
      }
    }

#if DEBUG_TMP
    std::cout<<"Finish fill position pairs.\n";
#endif

    std::partial_sort(matches_counts.begin(),
                      matches_counts.begin() + matches_counts.size() / 2,
                      matches_counts.end());
    size_t median_count = matches_counts[matches_counts.size() / 2 - 1];

    for (size_t i = 0; i < position_pairs.size(); i++)
    {
      auto itr_image_pair =
        matches.find(std::make_pair(epipolar_edges[i].first_id,
                                    epipolar_edges[i].second_id));
      if (itr_image_pair != matches.end())
      {
        position_pairs[i].weight =
          Scalar(itr_image_pair->second.size()) / Scalar(median_count);
      }
      else
      {
        position_pairs[i].weight = Scalar(1);
      }
    }

#if DEBUG_TMP
    std::cout<<"Finish compute position weight.\n";
#endif

    RotationContainer global_rotations;
    for (size_t image_id = 0; image_id < image_extrinsic_map.Size(); image_id++)
    {
      if (image_extrinsic_map.IsValid(image_id))
      {
        size_t extrinsic_id = image_extrinsic_map[image_id];
        global_rotations[image_id] =
          extrinsic_params_set[extrinsic_id].rotation();
      }
    }

#if DEBUG_TMP
    std::cout<<"Finish fill global rotation.\n";
#endif

    PositionAverage position_average;
    IndexContainer global_position_indices;
    PositionContainer global_positions;
    Err result = position_average(keysets,
                                  image_intrinsic_map,
                                  intrinsic_params_set,
                                  matches,
                                  position_pairs,
                                  global_rotations,
                                  global_position_indices,
                                  global_positions);

#if DEBUG_TMP
    std::cout<<"Finish position average.\n";
#endif

    if (result != 0) return result;

#if DEBUG_TMP
    std::cout<<"Copy global position.\n";
#endif

    hs::sfm::ObjectIndexMap image_extrinsic_map_position(
      image_extrinsic_map.Size());
    ExtrinsicParamsContainer extrinsic_params_set_position;
    for (size_t i = 0; i < global_position_indices.size(); i++)
    {
      size_t image_id = global_position_indices[i];
      image_extrinsic_map_position[image_id] = i;
      ExtrinsicParams extrinsic_params;
      size_t extrinsic_id = image_extrinsic_map[image_id];
      extrinsic_params.rotation() =
        extrinsic_params_set[extrinsic_id].rotation();
      extrinsic_params.position() = global_positions[i];
      extrinsic_params_set_position.push_back(extrinsic_params);
    }
    std::swap(image_extrinsic_map_position, image_extrinsic_map);
    std::swap(extrinsic_params_set_position, extrinsic_params_set);

#if DEBUG_TMP
    std::cout<<"TranslateAverage End\n";
#endif

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
