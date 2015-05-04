#ifndef _HS_SFM_SFM_PIPELINE_BUNDLE_ADJUSTMENT_OPTIMIZOR_HPP_
#define _HS_SFM_SFM_PIPELINE_BUNDLE_ADJUSTMENT_OPTIMIZOR_HPP_

#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_ceres_optimizor.hpp"

#define TMP_DEBUG 1
#if TMP_DEBUG
#include "hs_sfm/sfm_pipeline/reprojective_error_calculator.hpp"
#endif

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
class BundleAdjustmentOptimizor
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
  typedef hs::sfm::ba::CameraSharedVectorFunction<Scalar> BAVectorFunction;
  typedef typename BAVectorFunction::Index Index;
  typedef typename BAVectorFunction::XVector XVector;
  typedef typename BAVectorFunction::YVector YVector;
  typedef typename BAVectorFunction::FeatureMap FeatureMap;
  typedef typename BAVectorFunction::FeatureMapContainer FeatureMapContainer;
  typedef typename BAVectorFunction::ImageCameraMap ImageCameraMap;
  typedef typename BAVectorFunction::Point BAPoint;
  typedef typename BAVectorFunction::Image BAImage;
  typedef typename BAVectorFunction::Camera BACamera;

  typedef hs::sfm::ba::CameraSharedCeresOptimizor<BAVectorFunction>
          BAOptimizor;
  typedef typename BAOptimizor::YCovarianceInverse YCovarianceInverse;

  typedef typename IntrinsicParams::KMatrix KMatrix;
  typedef typename ExtrinsicParams::Position Position;

#if TMP_DEBUG
  typedef ReprojectiveErrorCalculator<Scalar> ReprojectiveErrorCalculatorType;
#endif

public:
  BundleAdjustmentOptimizor(size_t number_of_threads = 1)
    : number_of_threads_(number_of_threads) {}

  Err operator() (const ImageKeysetContainer& image_keysets,
                  const ImageIntrinsicMap& image_intrinsic_map,
                  const TrackContainer& tracks,
                  const ImageExtrinsicMap& image_extrinsic_map,
                  const TrackPointMap& track_point_map,
                  const ViewInfoIndexer& view_info_indexer,
                  IntrinsicParamsContainer& intrinsic_params_set,
                  ExtrinsicParamsContainer& extrinsic_params_set,
                  PointContainer& points) const
  {
#if TMP_DEBUG
    ReprojectiveErrorCalculatorType reprojective_error_calculator;
    Scalar error_before = reprojective_error_calculator(image_keysets,
                                                     intrinsic_params_set,
                                                     tracks,
                                                     image_intrinsic_map,
                                                     image_extrinsic_map,
                                                     track_point_map,
                                                     view_info_indexer,
                                                     extrinsic_params_set,
                                                     points);
    std::cout<<"error_before:"<<error_before<<"\n";
#endif
    size_t number_of_tracks = tracks.size();
    TrackContainer tracks_bundle;
    ObjectIndexMap bundle_point_map;
    Index number_of_keys = 0;
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
          number_of_keys += Index(track_bundle.size());
        }
      }
    }

    //设置feature_maps
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

    //设置image_camera_map
    ImageCameraMap image_camera_map;
    ObjectIndexMap extrinsic_image_map(extrinsic_params_set.size());
    for (size_t i = 0; i < image_extrinsic_map.Size(); i++)
    {
      if (image_extrinsic_map.IsValid(i))
      {
        extrinsic_image_map[image_extrinsic_map[i]] = i;
      }
    }
    for (size_t i = 0; i < extrinsic_image_map.Size(); i++)
    {
      if (!extrinsic_image_map.IsValid(i)) return -1;
      size_t image_id = extrinsic_image_map[i];
      image_camera_map.push_back(image_intrinsic_map[image_id]);
    }

    Index number_of_points = Index(number_of_tracks_bundle);
    Index number_of_images = Index(extrinsic_params_set.size());
    Index number_of_cameras = Index(intrinsic_params_set.size());
    BAVectorFunction vector_function;
    vector_function.set_number_of_images(number_of_images);
    vector_function.set_number_of_points(number_of_points);
    vector_function.set_number_of_keys(number_of_keys);
    vector_function.set_number_of_cameras(number_of_cameras);
    vector_function.set_feature_maps(feature_maps);
    vector_function.set_image_camera_map(image_camera_map);
    vector_function.intrinsic_computations_mask().set();

    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    XVector initial_x(x_size);
    for (Index i = 0; i < number_of_points; i++)
    {
      size_t point_id = bundle_point_map[size_t(i)];
      initial_x.template segment<BAVectorFunction::params_per_point_>(
        i * BAVectorFunction::params_per_point_) = points[point_id];
    }
    Index extrinsic_begin = vector_function.GetPointParamsSize();
    for (Index i = 0; i < number_of_images; i++)
    {
      const ExtrinsicParams& extrinsic_params = extrinsic_params_set[i];
      Point t = -(extrinsic_params.rotation() * extrinsic_params.position());
      Index offset = extrinsic_begin +
                     i * BAVectorFunction::extrinsic_params_per_image_;
      initial_x[offset + 0] = extrinsic_params.rotation()[0];
      initial_x[offset + 1] = extrinsic_params.rotation()[1];
      initial_x[offset + 2] = extrinsic_params.rotation()[2];
      initial_x[offset + 3] = t[0];
      initial_x[offset + 4] = t[1];
      initial_x[offset + 5] = t[2];
    }
    Index intrinsic_begin = vector_function.GetPointParamsSize() +
                            vector_function.GetExtrinsicParamsSize();
    Index intrinsic_params_size_per_camera =
      vector_function.GetIntrinsicParamsSizePerCamera();
    for (Index i = 0; i < number_of_cameras; i++)
    {
      Index offset = intrinsic_begin + i * intrinsic_params_size_per_camera;
      const IntrinsicParams& intrinsic_params = intrinsic_params_set[i];
      initial_x[offset + 0] = intrinsic_params.k1();
      initial_x[offset + 1] = intrinsic_params.k2();
      initial_x[offset + 2] = intrinsic_params.k3();
      initial_x[offset + 3] = intrinsic_params.d1();
      initial_x[offset + 4] = intrinsic_params.d2();
      initial_x[offset + 5] = intrinsic_params.focal_length();
      initial_x[offset + 6] = intrinsic_params.skew();
      initial_x[offset + 7] = intrinsic_params.principal_point_x();
      initial_x[offset + 8] = intrinsic_params.principal_point_y();
      initial_x[offset + 9] = intrinsic_params.pixel_ratio();
    }
    YVector near_y(y_size);
    Index feature_id = 0;
    for (size_t i = 0; i < number_of_tracks_bundle; i++)
    {
      size_t number_of_views_bundle = tracks_bundle[i].size();
      for (size_t j = 0; j < number_of_views_bundle; j++)
      {
        size_t image_id = tracks_bundle[i][j].first;
        size_t key_id = tracks_bundle[i][j].second;

        near_y.segment(feature_id * BAVectorFunction::params_per_key_,
                       BAVectorFunction::params_per_key_) =
          image_keysets[image_id][key_id];

        feature_id++;
      }
    }
    YCovarianceInverse y_covariance_inverse;
    y_covariance_inverse.SetKeysUniformStdDev(Scalar(1), number_of_keys);

    int max_num_iterations = 50;
    double function_tolerance = 1e-6;
    double parameter_tolerance = 1e-8;
    if (number_of_images > 2)
    {
      max_num_iterations = 10;
      function_tolerance = 1e-4;
      parameter_tolerance = 1e-4;
    }
#if DEBUG_TMP
    std::cout<<"Start optimizing.\n";
#endif
    BAOptimizor ba_optimizor(initial_x, int(number_of_threads_),
                             max_num_iterations, function_tolerance,
                             parameter_tolerance);
    XVector optimized_x;
    if (ba_optimizor(vector_function, near_y, y_covariance_inverse,
                     optimized_x) != 0)
    {
      return -1;
    }
#if DEBUG_TMP
    std::cout<<"Finish optimizing.\n";
#endif

    for (Index i = 0; i < number_of_points; i++)
    {
      size_t point_id = bundle_point_map[size_t(i)];
      points[point_id] = vector_function.GetPoint(i, optimized_x);
    }

    for (Index i = 0; i < number_of_images; i++)
    {
      ExtrinsicParams& extrinsic_params = extrinsic_params_set[i];
      BAImage ba_image = vector_function.GetImage(i, optimized_x);
      extrinsic_params.rotation()[0] = ba_image[0];
      extrinsic_params.rotation()[1] = ba_image[1];
      extrinsic_params.rotation()[2] = ba_image[2];
      Position position = ba_image.template segment<3>(3);
      extrinsic_params.position() =
        -(extrinsic_params.rotation().Inverse() * position);
    }

    for (Index i = 0; i < number_of_cameras; i++)
    {
      IntrinsicParams& intrinsic_params = intrinsic_params_set[i];
      BACamera ba_camera = vector_function.GetCamera(i, optimized_x);
      intrinsic_params.set_k1(ba_camera[0]);
      intrinsic_params.set_k2(ba_camera[1]);
      intrinsic_params.set_k3(ba_camera[2]);
      intrinsic_params.set_d1(ba_camera[3]);
      intrinsic_params.set_d2(ba_camera[4]);
      intrinsic_params.set_focal_length(ba_camera[5]);
      intrinsic_params.set_skew(ba_camera[6]);
      intrinsic_params.set_principal_point_x(ba_camera[7]);
      intrinsic_params.set_principal_point_y(ba_camera[8]);
      intrinsic_params.set_pixel_ratio(ba_camera[9]);
#if TMP_DEBUG
      std::cout<<"camera "<<i<<":\n";
      std::cout<<"k1:"<<ba_camera[0]<<"\n";
      std::cout<<"k2:"<<ba_camera[1]<<"\n";
      std::cout<<"k3:"<<ba_camera[2]<<"\n";
      std::cout<<"d1:"<<ba_camera[3]<<"\n";
      std::cout<<"d2:"<<ba_camera[4]<<"\n";
      std::cout<<"focal_length:"<<ba_camera[5]<<"\n";
      std::cout<<"skew:"<<ba_camera[6]<<"\n";
      std::cout<<"principal_x:"<<ba_camera[7]<<"\n";
      std::cout<<"principal_y:"<<ba_camera[8]<<"\n";
      std::cout<<"pixel_ratio:"<<ba_camera[9]<<"\n";
#endif
    }

#if TMP_DEBUG
    Scalar error_after = reprojective_error_calculator(image_keysets,
                                                       intrinsic_params_set,
                                                       tracks,
                                                       image_intrinsic_map,
                                                       image_extrinsic_map,
                                                       track_point_map,
                                                       view_info_indexer,
                                                       extrinsic_params_set,
                                                       points);
    std::cout<<"error_after:"<<error_after<<"\n";
    std::cout<<"number_of_points:"<<number_of_points<<"\n";
#endif
    return 0;
  }

private:
  size_t number_of_threads_;
};

}
}
}

#endif
