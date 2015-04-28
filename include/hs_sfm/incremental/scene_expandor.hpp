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

#define DEBUG_TMP 1
#if DEBUG_TMP
#include <sstream>
#include <fstream>
#include <iomanip>
#include "hs_sfm/incremental/reprojective_error_calculator.hpp"
#include "hs_sfm/sfm_file_io/scene_ply_saver.hpp"
#include "hs_sfm/sfm_utility/debug_tmp.hpp"
#include "hs_sfm/sfm_utility/similar_transform_estimator.hpp"
#endif

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

#if DEBUG_TMP
  typedef DebugTrue<Scalar> DebugTrueType;
#endif

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
                  hs::progress::ProgressManager* progress_manager = NULL
#if DEBUG_TMP
                  , const DebugTrueType& debug_true = DebugTrueType()
#endif
                  ) const
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
               progress_manager
#if DEBUG_TMP
               , debug_true
#endif
               );
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
          hs::progress::ProgressManager* progress_manager = NULL
#if DEBUG_TMP
          , const DebugTrueType& debug_true = DebugTrueType()
#endif
          ) const
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

#if DEBUG_TMP
    double total_bundle_adjustment_time = 0.0;
    double total_image_expansion_time = 0.0;
    double total_point_expansion_time = 0.0;

    std::string points_added_accuracy_path =
      debug_true.prefix() + "_points_added_accuracy.txt";
    std::ofstream points_added_accuracy_file(points_added_accuracy_path);
#endif

    while (1)
    {
      if (progress_manager)
      {
        if (!progress_manager->CheckKeepWorking())
        {
          break;
        }
      }

#if DEBUG_TMP
      ReprojectiveErrorCalculator<Scalar> reprojective_error_calculator;
      Scalar reprojective_error_before =
        reprojective_error_calculator(image_keysets,
                                      intrinsic_params_set,
                                      tracks,
                                      image_intrinsic_map,
                                      image_extrinsic_map,
                                      track_point_map,
                                      view_info_indexer,
                                      extrinsic_params_set,
                                      points);

      std::chrono::time_point<std::chrono::system_clock> start, end;
      start = std::chrono::system_clock::now();

      std::stringstream ss;
      ss<<"extrinsic_before_"<<extrinsic_params_set.size()<<".ply";
      std::string extrinsic_before_path;
      ss>>extrinsic_before_path;
      OutputScene(extrinsic_before_path,
                  extrinsic_params_set, points, intrinsic_params_set[0]);
#endif
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
#if DEBUG_TMP
      end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      total_bundle_adjustment_time += elapsed_seconds.count();
      std::cout<<"Bundle Adjustment took "<<elapsed_seconds.count()
               <<" seconds.\n";
      Scalar reprojective_error_after =
        reprojective_error_calculator(image_keysets,
                                      intrinsic_params_set,
                                      tracks,
                                      image_intrinsic_map,
                                      image_extrinsic_map,
                                      track_point_map,
                                      view_info_indexer,
                                      extrinsic_params_set,
                                      points);
      ss.clear();
      ss.str("");
      ss<<"intrinsic_"<<extrinsic_params_set.size()<<".txt";
      std::string intrinsic_path;
      ss>>intrinsic_path;
      std::ofstream intrinsic_file(intrinsic_path.c_str());
      intrinsic_file.setf(std::ios::fixed);
      intrinsic_file<<std::setprecision(8);

      intrinsic_file<<points.size()<<"\n";
      intrinsic_file<<reprojective_error_before<<" "
                    <<reprojective_error_after<<"\n";
      for (size_t i = 0; i < intrinsic_params_set.size(); i++)
      {
        const IntrinsicParams& intrinsic_param = intrinsic_params_set[i];
        intrinsic_file << i << " "
                       << intrinsic_param.focal_length() << " "
                       << intrinsic_param.skew() << " "
                       << intrinsic_param.principal_point_x() << " "
                       << intrinsic_param.principal_point_y() << " "
                       << intrinsic_param.pixel_ratio() << " "
                       << intrinsic_param.k1() << " "
                       << intrinsic_param.k2() << " "
                       << intrinsic_param.k3() << " "
                       << intrinsic_param.d1() << " "
                       << intrinsic_param.d2() << "\n";
      }

      intrinsic_file.close();

      ss.clear();
      ss.str("");
      ss<<"extrinsic_after_"<<extrinsic_params_set.size()<<".ply";
      std::string extrinsic_after_path;
      ss>>extrinsic_after_path;
      OutputScene(extrinsic_after_path,
                  extrinsic_params_set, points, intrinsic_params_set[0]);

#endif

#if DEBUG_TMP
      start = std::chrono::system_clock::now();
#endif
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
#if DEBUG_TMP
        std::cout<<"Added Image "<<new_image_ids[i]<<"\n";
        std::cout<<"Position:\n"<<new_extrinsic_params_set[i].position()<<"\n";
#endif
        extrinsic_params_set.push_back(new_extrinsic_params_set[i]);
        image_extrinsic_map[new_image_ids[i]] = extrinsic_params_set.size() - 1;
      }

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

#if DEBUG_TMP
      std::cout<<new_image_ids.size()
               <<" images added(Total: "<<extrinsic_params_set.size()<<").\n";
      end = std::chrono::system_clock::now();
      elapsed_seconds = end - start;
      total_image_expansion_time += elapsed_seconds.count();

      start = std::chrono::system_clock::now();
      size_t number_of_points_before_adding = points.size();
#endif

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

#if DEBUG_TMP
      end = std::chrono::system_clock::now();
      elapsed_seconds = end - start;
      total_point_expansion_time += elapsed_seconds.count();

      if (!debug_true.points_true().empty())
      {
        PointContainer points_added_estimate;
        PointContainer points_added_true;
        for (size_t i = 0; i < tracks.size(); i++)
        {
          if (track_point_map.IsValid(i))
          {
            size_t point_id = track_point_map[i];
            if (point_id >= number_of_points_before_adding)
            {
              points_added_estimate.push_back(points[point_id]);
              points_added_true.push_back(debug_true.points_true()[i]);
            }
          }
        }

        SimilarTransformEstimator<Scalar> similar_estimator;
        typename SimilarTransformEstimator<Scalar>::Rotation rotation_similar;
        typename SimilarTransformEstimator<Scalar>::Translate translate_similar;
        Scalar scale_similar;
        similar_estimator(points_added_estimate, points_added_true,
                          rotation_similar, translate_similar, scale_similar);
        Scalar error_planar = Scalar(0);
        Scalar error_height = Scalar(0);
        for (size_t i = 0; i < points_added_estimate.size(); i++)
        {
          Point point_estimate = points_added_estimate[i];
          const Point& point_true = points_added_true[i];
          point_estimate = scale_similar * (rotation_similar * point_estimate) +
                           translate_similar;
          error_planar += (point_estimate.template segment<2>(0) -
                           point_true.template segment<2>(0)).norm();
          error_height += std::abs(point_estimate[2] - point_true[2]);
        }
        error_planar /= Scalar(points_added_estimate.size());
        error_height /= Scalar(points_added_estimate.size());

        points_added_accuracy_file<<extrinsic_params_set.size()<<" "
                                  <<points_added_estimate.size()<<" "
                                  <<error_planar<<" "
                                  <<error_height<<"\n";
      }

#endif

      if (progress_manager)
      {
        progress_manager->SetCurrentSubProgressCompleteRatio(
          float(extrinsic_params_set.size() / float(image_keysets.size())));
      }
    }

#if DEBUG_TMP
    std::ofstream time_file("time.txt");
    time_file<<"Bundle adjustment took "<<total_bundle_adjustment_time<<" seconds\n";
    time_file<<"Image expansion took "<<total_image_expansion_time<<" seconds\n";
    time_file<<"Point expansion took "<<total_point_expansion_time<<" seconds\n";
    points_added_accuracy_file.close();
#endif

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

#if DEBUG_TMP
  void OutputScene(const std::string& scene_path,
                   const ExtrinsicParamsContainer& extrinsic_params_set,
                   const PointContainer& points,
                   const IntrinsicParams& intrinsic_params) const
  {
    typedef hs::sfm::fileio::ScenePLYSaver<Scalar, size_t> SceneSaver;
    typedef typename SceneSaver::Image Image;
    typedef typename SceneSaver::ImageContainer ImageContainer;

    IntrinsicParamsContainer intrinsic_params_out(extrinsic_params_set.size(),
                                                  intrinsic_params);
    Image image_out;
    image_out.m_width = 6000;
    image_out.m_height = 4000;
    ImageContainer images_out(extrinsic_params_set.size(), image_out);
    SceneSaver saver(Scalar(0.5));
    saver(scene_path,
          intrinsic_params_out, extrinsic_params_set, images_out, points);
  }
#endif

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
