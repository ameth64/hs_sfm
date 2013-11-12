#ifndef _HS_SFM_INCREMENTAL_SCENE_EXPANDOR_HPP_
#define _HS_SFM_INCREMENTAL_SCENE_EXPANDOR_HPP_

#include <vector>
#include <map>

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/projective/single_camera_params_maximum_likelihood_estimator.hpp"
#include "hs_sfm/projective/pmatrix_dlt_ransac_refiner.hpp"
#include "hs_sfm/triangulate/multiple_view_maximum_likelihood_estimator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_levenberg_marquardt_optimizor.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

class ObjectIndexMap
{
public:
  ObjectIndexMap(){}
  ObjectIndexMap(size_t number_of_objects)
    : mapper_(number_of_objects, invalid_value()) {}
public:

  size_t GetMappedId(size_t object_id) const
  {
    return mapper_[object_id];
  }

  void SetObjectId(size_t object_id, size_t mapped_id)
  {
    mapper_[object_id] = mapped_id;
  }

  size_t operator[] (size_t object_id) const
  {
    return mapper_[object_id];
  }

  size_t& operator[] (size_t object_id)
  {
    return mapper_[object_id];
  }

  bool IsValid(size_t object_id) const
  {
    return (mapper_[object_id] != invalid_value());
  }

  void Resize(size_t number_of_objects)
  {
    mapper_.resize(number_of_objects, invalid_value());
  }

  size_t Size() const
  {
    return mapper_.size();
  }

  void AddObject(size_t mapped_id)
  {
    mapper_.push_back(mapped_id);
  }

private:
  inline size_t invalid_value() const
  {
    return std::numeric_limits<size_t>::max();
  }

private:
  std::vector<size_t> mapper_;
};

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
  struct ViewInfo
  {
    size_t track_id;
    size_t image_id;
    size_t key_id;
    bool is_blunder;
  };

  class ViewInfoIndexer
  {
  public:
    void SetViewInfoByTracks(const hs::sfm::TrackContainer& tracks,
                             const TrackPointMap& track_point_map,
                             const ImageExtrinsicMap& image_extrinsic_map)
    {
      views_info_.clear();
      track_image_index_.clear();
      image_key_index_.clear();
      size_t number_of_tracks = tracks.size();
      for (size_t i = 0; i < number_of_tracks; i++)
      {
        size_t number_of_views = tracks[i].size();
        for (size_t j = 0; j < number_of_views; j++)
        {
          ViewInfo view_info;
          view_info.track_id = i;
          view_info.image_id = tracks[i][j].first;
          view_info.key_id = tracks[i][j].second;
          view_info.is_blunder = false;
          views_info_.push_back(view_info);
          size_t view_id = views_info_.size() - 1;
          track_image_index_[std::make_pair(view_info.track_id,
                                            view_info.image_id)] = view_id;
          image_key_index_[std::make_pair(view_info.image_id,
                                          view_info.key_id)] = view_id;
        }
      }
    }

    const ViewInfo& GetViewInfoByTrackImage(size_t track_id,
                                            size_t image_id) const
    {
      auto view_itr = track_image_index_.find(std::make_pair(track_id,
                                                              image_id));
      return views_info_[view_itr->second];
    }

    ViewInfo& GetViewInfoByTrackImage(size_t track_id,
                                      size_t image_id)
    {
      return views_info_[track_image_index_[std::make_pair(track_id,
                                                           image_id)]];
    }

    const ViewInfo& GetViewInfoByImageKey(size_t image_id,
                                        size_t key_id) const
    {
      return views_info_[image_key_index_[std::make_pair(image_id, key_id)]];
    }

    ViewInfo& GetViewInfoByImageKey(size_t image_id,
                                  size_t key_id)
    {
      return views_info_[image_key_index_[std::make_pair(image_id, key_id)]];
    }

  private:
    std::vector<ViewInfo> views_info_;
    std::map<std::pair<size_t, size_t>, size_t> track_image_index_;
    std::map<std::pair<size_t, size_t>, size_t> image_key_index_;
  };

  typedef std::vector<size_t> ImageViewTracks;
  typedef std::vector<ImageViewTracks> ImageViewTracksContainer;

  typedef
    hs::sfm::projective::SingleCameraParamsMaximumLikelihoodEstimator<Scalar>
      CameraParamsMLEstimator;
  typedef typename CameraParamsMLEstimator::Correspondence Correspondence;
  typedef typename CameraParamsMLEstimator::CorrespondenceContainer
                   CorrespondenceContainer;
  typedef typename CameraParamsMLEstimator::KeyCovariance KeyCovariance;

  typedef hs::sfm::projective::PMatrixDLTRansacRefiner<Scalar>
          PMatrixRasacRefiner;
  typedef typename PMatrixRasacRefiner::IndexSet IndexSet;

  typedef hs::sfm::triangulate::MultipleViewMaximumLikelihoodEstimator<Scalar>
    TriangulateMLEstimator;
  typedef typename TriangulateMLEstimator::KeyContainer KeyContainer;

  typedef CameraFunctions<Scalar> CamFunctions;
  typedef typename CamFunctions::ProjectionMatrix PMatrix;

  typedef hs::sfm::ba::BANaiveVectorFunction<Scalar> BAVectorFunction;
  typedef typename BAVectorFunction::Index Index;
  typedef typename BAVectorFunction::XVector XVector;
  typedef typename BAVectorFunction::YVector YVector;
  typedef typename BAVectorFunction::FeatureMap FeatureMap;
  typedef typename BAVectorFunction::FeatureMapContainer FeatureMapContainer;

  typedef hs::sfm::ba::BANaiveLevenbergMarquardtOptimizor<BAVectorFunction>
          BAOptimizor;
  typedef typename BAOptimizor::YCovarianceInverse YCovarianceInverse;

  typedef typename IntrinsicParams::KMatrix KMatrix;

public:
  SceneExpandor(
    size_t add_new_image_matches_threshold = 8,
    Scalar pmatrix_ransac_threshold = 4.0,
    size_t min_triangulate_views = 2,
    Scalar triangulate_error_threshold = 4.0)
    : add_new_image_matches_threshold_(add_new_image_matches_threshold),
      pmatrix_ransac_threshold_(pmatrix_ransac_threshold),
      min_triangulate_views_(min_triangulate_views),
      triangulate_error_threshold_(triangulate_error_threshold){}

  Err operator() (const ImageKeysetContainer& image_keysets,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const TrackContainer& tracks,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  ImageExtrinsicMap& image_extrinsic_map,
                  PointContainer& points,
                  TrackPointMap& track_point_map,
                  Scalar& reprojection_error) const
  {
    return Run(image_keysets,
               intrinsic_params_set,
               tracks,
               extrinsic_params_set,
               image_extrinsic_map,
               points,
               track_point_map,
               reprojection_error);
  }

  Err Run(const ImageKeysetContainer& image_keysets,
          const IntrinsicParamsContainer& intrinsic_params_set,
          const TrackContainer& tracks,
          ExtrinsicParamsContainer& extrinsic_params_set,
          ImageExtrinsicMap& image_extrinsic_map,
          PointContainer& points,
          TrackPointMap& track_point_map,
          Scalar& reprojection_error) const
  {
    ImageViewTracksContainer image_view_tracks_set;
    ViewInfoIndexer view_info_indexer;
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
      Scalar reprojection_error_before =
      GetReprojectiveError(image_keysets,
                           intrinsic_params_set,
                           tracks,
                           image_extrinsic_map,
                           track_point_map,
                           view_info_indexer,
                           extrinsic_params_set,
                           points);
      reprojection_error = reprojection_error_before;
      if (BundleAdjustment(image_keysets,
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
      reprojection_error =
        GetReprojectiveError(image_keysets,
                             intrinsic_params_set,
                             tracks,
                             image_extrinsic_map,
                             track_point_map,
                             view_info_indexer,
                             extrinsic_params_set,
                             points);

      ExtrinsicParams new_extrinsic_params;
      size_t new_image_id;
      if (AddNewImage(image_keysets,
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

      if (TriangulateNewPoints(image_keysets,
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
    view_info_indexer.SetViewInfoByTracks(tracks,
                                          track_point_map,
                                          image_extrinsic_map);

    return 0;
  }

  Err AddNewImage(const ImageKeysetContainer& image_keysets,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const PointContainer& points,
                  const TrackPointMap& track_point_map,
                  const ImageViewTracksContainer& image_view_tracks_set,
                  const ImageExtrinsicMap& image_extrinsic_map,
                  ViewInfoIndexer& view_info_indexer,
                  ExtrinsicParams& new_extrinsic_params,
                  size_t& new_image_id) const
  {
    size_t number_of_images = image_keysets.size();
    if (number_of_images != intrinsic_params_set.size())
    {
      return -1;
    }

    //查找拍摄到当前点云的其他影像，选取匹配数最多的影像，
    //若匹配数最多的影像的匹配数没超过阈值，则函数执行不成功。
    size_t max_matches_image_id = 0;
    size_t max_number_of_view_tracks = 0;
    ImageViewTracks max_image_view_tracks;
    for (size_t i = 0; i < number_of_images; i++)
    {
      //已加入的影像排除在外
      if (image_extrinsic_map.IsValid(i))
      {
        continue;
      }

      size_t number_of_view_tracks = 0;
      ImageViewTracks image_view_tracks;
      for (size_t j = 0; j < image_view_tracks_set[i].size(); j++)
      {
        if (!track_point_map.IsValid(image_view_tracks_set[i][j]))
        {
          continue;
        }

        image_view_tracks.push_back(image_view_tracks_set[i][j]);
        number_of_view_tracks++;
      }

      if (number_of_view_tracks > max_number_of_view_tracks)
      {
        max_number_of_view_tracks = number_of_view_tracks;
        max_matches_image_id = i;
        max_image_view_tracks.swap(image_view_tracks);
      }
    }

    if (max_number_of_view_tracks < add_new_image_matches_threshold_)
    {
      return -1;
    }

    new_image_id = max_matches_image_id;

    //获取三维点与二维点的对应
    CorrespondenceContainer coarse_correspondences;
    ObjectIndexMap key_track_map(max_number_of_view_tracks);
    for (size_t i = 0; i < max_number_of_view_tracks; i++)
    {
      size_t track_id = max_image_view_tracks[i];
      if (!track_point_map.IsValid(track_id))
      {
        return -1;
      }
      size_t point_id = track_point_map[track_id];

      ViewInfo& view_info =
        view_info_indexer.GetViewInfoByTrackImage(track_id, new_image_id);
      view_info.is_blunder = true;

      Correspondence correspondence;
      correspondence.first = image_keysets[new_image_id][view_info.key_id];
      correspondence.second = points[point_id];
      coarse_correspondences.push_back(correspondence);
      key_track_map[i] = track_id;
    }

    //Ransac剔除错误的对应
    IndexSet inlier_indices;
    PMatrixRasacRefiner ransac_refiner;
    CorrespondenceContainer refined_correspondences;
    if (ransac_refiner(coarse_correspondences,
                       pmatrix_ransac_threshold_,
                       refined_correspondences,
                       inlier_indices) != 0)
    {
      return -1;
    }
    size_t number_of_inliers = inlier_indices.size();
    for (size_t i = 0; i < number_of_inliers; i++)
    {
      size_t track_id = key_track_map[inlier_indices[i]];
      ViewInfo& view_info =
        view_info_indexer.GetViewInfoByTrackImage(track_id, new_image_id);
      view_info.is_blunder = false;
    }

    //使用剔除粗差的点计算外方位元素
    CameraParamsMLEstimator camera_ml_estimator;
    KeyCovariance key_covariance = KeyCovariance::Identity();
    IntrinsicParams intrinsic_params_estimate;
    if (camera_ml_estimator(refined_correspondences,
                            key_covariance,
                            intrinsic_params_set[new_image_id],
                            Scalar(1),
                            Scalar(0.01),
                            Scalar(1),
                            Scalar(1),
                            Scalar(1e-4),
                            intrinsic_params_estimate,
                            new_extrinsic_params) != 0)
    {
      return -1;
    }

    return 0;
  }

  Err TriangulateNewPoints(const ImageKeysetContainer& image_keysets,
                           const IntrinsicParamsContainer& intrinsic_params_set,
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
            size_t extrinsic_id = image_extrinsic_map[image_id];
            intrinsic_params_set_view[i] =
              intrinsic_params_set[image_id];
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
              PMatrix P =
                CamFunctions::GetProjectionMatrix(intrinsic_params_set_view[i],
                                                  extrinsic_params_set_view[i]);
              EIGEN_VECTOR(Scalar, 3) hkey = P.block(0, 0, 3, 3) * point +
                                             P.col(3);
              hkey /= hkey(2);
              Scalar error = (keys[i] - hkey.segment(0, 2)).norm();

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

  Err BundleAdjustment(const ImageKeysetContainer& image_keysets,
                       const IntrinsicParamsContainer& intrinsic_params_set,
                       const TrackContainer& tracks,
                       const ImageExtrinsicMap& image_extrinsic_map,
                       const TrackPointMap& track_point_map,
                       const ViewInfoIndexer& view_info_indexer,
                       ExtrinsicParamsContainer& extrinsic_params_set,
                       PointContainer& points) const
  {
    size_t number_of_tracks = tracks.size();
    TrackContainer tracks_bundle;
    ObjectIndexMap bundle_point_map;
    Index number_of_features = 0;
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
            track_bundle.push_back(tracks[track_id][i]);
          }
        }
        if (!track_bundle.empty())
        {
          tracks_bundle.push_back(track_bundle);
          size_t point_id = track_point_map[track_id];
          bundle_point_map.AddObject(point_id);
          number_of_features += Index(track_bundle.size());
        }
      }
    }

    size_t number_of_tracks_bundle = tracks_bundle.size();
    FeatureMapContainer feature_maps;
    for (size_t i = 0; i < number_of_tracks_bundle; i++)
    {
      size_t number_of_views_bundle = tracks_bundle[i].size();
      for (size_t j = 0; j < number_of_views_bundle; j++)
      {
        size_t image_id = tracks_bundle[i][j].first;
        size_t key_id = tracks_bundle[i][j].second;
        size_t extrinsic_id = image_extrinsic_map[image_id];

        feature_maps.push_back(std::make_pair(extrinsic_id, i));
      }
    }

    Index number_of_points = Index(number_of_tracks_bundle);
    Index number_of_cameras = Index(extrinsic_params_set.size());
    BAVectorFunction vector_function;
    vector_function.set_number_of_cameras(number_of_cameras);
    vector_function.set_number_of_points(number_of_points);
    vector_function.set_number_of_features(number_of_features);
    vector_function.set_feature_maps(feature_maps);

    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    XVector initial_x(x_size);
    YVector near_y(y_size);
    for (Index i = 0; i < number_of_cameras; i++)
    {
      const ExtrinsicParams& extrinsic_params = extrinsic_params_set[i];
      Point t = -(extrinsic_params.rotation() * extrinsic_params.position());
      initial_x[i * BAVectorFunction::params_per_camera_ + 0] =
        extrinsic_params.rotation()[0];
      initial_x[i * BAVectorFunction::params_per_camera_ + 1] =
        extrinsic_params.rotation()[1];
      initial_x[i * BAVectorFunction::params_per_camera_ + 2] =
        extrinsic_params.rotation()[2];
      initial_x[i * BAVectorFunction::params_per_camera_ + 3] = t[0];
      initial_x[i * BAVectorFunction::params_per_camera_ + 4] = t[1];
      initial_x[i * BAVectorFunction::params_per_camera_ + 5] = t[2];
    }
    Index camera_params_size = vector_function.GetCameraParamsSize();
    for (Index i = 0; i < number_of_points; i++)
    {
      size_t point_id = bundle_point_map[size_t(i)];
      initial_x.segment(
        camera_params_size + i * BAVectorFunction::params_per_point_,
        BAVectorFunction::params_per_point_) = points[point_id];
    }
    Index feature_id = 0;
    for (size_t i = 0; i < number_of_tracks_bundle; i++)
    {
      size_t number_of_views_bundle = tracks_bundle[i].size();
      for (size_t j = 0; j < number_of_views_bundle; j++)
      {
        size_t image_id = tracks_bundle[i][j].first;
        size_t key_id = tracks_bundle[i][j].second;
        size_t extrinsic_id = image_extrinsic_map[image_id];

        KMatrix K = intrinsic_params_set[image_id].GetKMatrix();
        EIGEN_VECTOR(Scalar, 3) hkey;
        hkey.template segment<2>(0) = image_keysets[image_id][key_id];
        hkey(2) = Scalar(1);
        hkey = K.inverse() * hkey;
        hkey /= hkey(2);

        near_y.segment(feature_id * BAVectorFunction::params_per_feature_,
                       BAVectorFunction::params_per_feature_) =
          hkey.segment(0, 2);
        feature_id++;
      }
    }
    YCovarianceInverse y_covariance_inverse;
    for (Index i = 0; i < number_of_features; i++)
    {
      y_covariance_inverse.blocks.push_back(
        EIGEN_MATRIX(Scalar, 2, 2)::Identity());
    }

    BAOptimizor ba_optimizor(initial_x);
    XVector optimized_x;
    if (ba_optimizor(vector_function, near_y, y_covariance_inverse,
                     optimized_x) != 0)
    {
      return -1;
    }

    for (Index i = 0; i < number_of_cameras; i++)
    {
      ExtrinsicParams& extrinsic_params = extrinsic_params_set[i];
      extrinsic_params.rotation()[0] =
        optimized_x[i * BAVectorFunction::params_per_camera_ + 0];
      extrinsic_params.rotation()[1] =
        optimized_x[i * BAVectorFunction::params_per_camera_ + 1];
      extrinsic_params.rotation()[2] =
        optimized_x[i * BAVectorFunction::params_per_camera_ + 2];
      Point t =
        optimized_x.segment(i * BAVectorFunction::params_per_camera_ + 3, 3);
      Point c = -((extrinsic_params.rotation().Inverse()) * t);
      extrinsic_params.position() = c;
    }

    for (Index i = 0; i < number_of_points; i++)
    {
      size_t point_id = bundle_point_map[size_t(i)];
      points[point_id] =
        optimized_x.segment(
        camera_params_size + i * BAVectorFunction::params_per_point_,
        BAVectorFunction::params_per_point_);
    }
    
    return 0;
  }

  Scalar GetReprojectiveError(const ImageKeysetContainer& image_keysets,
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
            PMatrix P = CamFunctions::GetProjectionMatrix(intrinsic_params,
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

protected:
  size_t add_new_image_matches_threshold_;
  Scalar pmatrix_ransac_threshold_;
  size_t min_triangulate_views_;
  Scalar triangulate_error_threshold_;
};

}
}
}

#endif
