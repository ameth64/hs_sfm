#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_SYNTHETIC_DATA_GENERATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_SYNTHETIC_DATA_GENERATOR_HPP_

#include "boost/random.hpp"

#include "hs_sfm/synthetic/multiple_flight_generator.hpp"
#include "hs_sfm/synthetic/multiple_camera_keyset_generator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_common_types.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Integer>
class IntegerUniformRandomGenerator
{
public:
  typedef _Integer Integer;
  typedef int Err;

  static Err Generate(Integer min, Integer max, Integer& var)
  {
    if (max < min)
    {
      return -1;
    }

    boost::random::uniform_int_distribution<Integer> distribution(min, max);
    var = distribution(random_generator_);

    return 0;
  }

private:
  static boost::random::mt19937 random_generator_;
};
template <typename _Integer>
boost::random::mt19937
IntegerUniformRandomGenerator<_Integer>::random_generator_;


template <typename _Scalar, typename _ImageDimension>
class CameraSharedSyntheticDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

  typedef hs::sfm::synthetic::MultipleFlightGenerator<Scalar, ImageDimension>
          MultipleFlightGenerator;
  typedef typename MultipleFlightGenerator::FlightGenerator FlightGenerator;
  typedef typename MultipleFlightGenerator::FlightGeneratorContainer
                   FlightGeneratorContainer;
  typedef typename MultipleFlightGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename MultipleFlightGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename MultipleFlightGenerator::Image Image;
  typedef typename MultipleFlightGenerator::ImageContainer ImageContainer;
  typedef typename MultipleFlightGenerator::Point3D Point3D;
  typedef typename MultipleFlightGenerator::Point3DContainer Point3DContainer;

  typedef hs::sfm::synthetic::MultipleCameraKeysetGenerator<Scalar,
                                                            ImageDimension>
          KeysetGenerator;
  typedef typename KeysetGenerator::IntrinsicParams IntrinsicParams;
  typedef typename KeysetGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::FeatureMap FeatureMap;
  typedef typename VectorFunction::FeatureMapContainer FeatureMapContainer;
  typedef typename VectorFunction::ImageCameraMap ImageCameraMap;
  typedef typename VectorFunction::Point BAPoint;
  typedef typename VectorFunction::PointContainer BAPointContainer;
  typedef typename VectorFunction::Image BAImage;
  typedef typename VectorFunction::ImageContainer BAImageContainer;
  typedef typename VectorFunction::Camera BACamera;
  typedef typename VectorFunction::CameraContainer BACameraContainer;

private:
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;

public:
  CameraSharedSyntheticDataGenerator(
    Scalar flight_longitudinal_overlap_ratio,
    Scalar flight_lateral_overlap_ratio,
    Scalar north_west_angle,
    Scalar north_west_angle_stddev,
    Scalar offset_stddev,
    const FlightGeneratorContainer& flight_generators,
    size_t number_of_points,
    size_t number_of_planar_constrained_points,
    size_t number_of_full_constrained_points,
    size_t number_of_constrained_images,
    size_t number_of_constrained_cameras,
    const IntrinsicParamsContainer& intrinsic_params_set,
    const FixMask& fix_mask)
    : multiple_flight_generator_(flight_longitudinal_overlap_ratio,
                                 flight_lateral_overlap_ratio,
                                 north_west_angle,
                                 north_west_angle_stddev,
                                 offset_stddev,
                                 flight_generators),
      keyset_generator_(),
      number_of_points_(number_of_points),
      number_of_planar_constrained_points_(number_of_planar_constrained_points),
      number_of_full_constrained_points_(number_of_full_constrained_points),
      number_of_constrained_images_(number_of_constrained_images),
      number_of_constrained_cameras_(number_of_constrained_cameras),
      intrinsic_params_set_(intrinsic_params_set),
      fix_mask_(fix_mask) {}

  Err operator() (VectorFunction& vector_function,
                  XVector& x,
                  YVector& y) const
  {
    Point3DContainer points;
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    ImageCameraMap image_camera_map;
    KeysetContainer keysets;
    TrackContainer tracks;
    CameraViewContainer camera_views;
    if (GenerateSyntheticData(points,
                              extrinsic_params_set,
                              images,
                              image_camera_map,
                              keysets,
                              tracks,
                              camera_views) != 0) return -1;

    std::vector<Index> keyset_id_offsets;
    std::vector<size_t> image_map;
    std::vector<size_t> point_map;
    size_t number_of_keys;
    size_t number_of_available_images;
    size_t number_of_available_points;
    FeatureMapContainer feature_maps;
    if (GenerateMaps(tracks,
                     camera_views,
                     keysets,
                     keyset_id_offsets,
                     image_map,
                     point_map,
                     number_of_keys,
                     number_of_available_images,
                     number_of_available_points,
                     feature_maps) != 0)
    {
      return -1;
    }

    std::vector<size_t> constrained_point_ids;
    std::vector<size_t> constrained_image_ids;
    if (GenerateVectorFunction(number_of_available_images,
                               number_of_available_points,
                               number_of_keys,
                               image_map,
                               point_map,
                               tracks,
                               camera_views,
                               feature_maps,
                               image_camera_map,
                               points,
                               extrinsic_params_set,
                               constrained_image_ids,
                               constrained_point_ids,
                               vector_function) != 0)
    {
      return -1;
    }

    if (GenerateXVector(points,
                        tracks,
                        extrinsic_params_set,
                        camera_views,
                        vector_function,
                        x) != 0)
    {
      return -1;
    }

    if (GenerateYVector(keysets,
                        extrinsic_params_set,
                        points,
                        tracks,
                        vector_function,
                        constrained_image_ids,
                        constrained_point_ids,
                        number_of_available_points,
                        y) != 0)
    {
      return -1;
    }

    return 0;
  }

private:
  Err GenerateSyntheticData(Point3DContainer& points,
                            ExtrinsicParamsContainer& extrinsic_params_set,
                            ImageContainer& images,
                            ImageCameraMap& image_camera_map,
                            KeysetContainer& keysets,
                            TrackContainer& tracks,
                            CameraViewContainer& camera_views) const
  {
    points.clear();
    if (multiple_flight_generator_.GeneratePoints(number_of_points_,
                                                  points) != 0) return -1;

    size_t number_of_cameras = intrinsic_params_set_.size();
    if (number_of_cameras != multiple_flight_generator_.GetNumberOfFlights())
    {
      return -1;
    }

    extrinsic_params_set.clear();
    images.clear();
    image_camera_map.clear();
    std::vector<size_t> image_intrinsic_map;
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      ExtrinsicParamsContainer extrinsic_params_set_flight;
      ImageContainer images_flight;
      if (multiple_flight_generator_.GenerateExtrinsicParamsContainer(
            i, extrinsic_params_set_flight, images_flight) != 0) return -1;

      for (size_t j = 0; j < images_flight.size(); j++)
      {
        extrinsic_params_set.push_back(extrinsic_params_set_flight[j]);
        images.push_back(images_flight[j]);
        image_camera_map.push_back(Index(i));
        image_intrinsic_map.push_back(i);
      }
    }

    keysets.clear();
    tracks.clear();
    camera_views.clear();
    if (keyset_generator_(intrinsic_params_set_,
                          extrinsic_params_set,
                          images,
                          points,
                          image_intrinsic_map,
                          keysets,
                          tracks,
                          camera_views) != 0) return -1;

    return 0;
  }

  Err GenerateMaps(const TrackContainer& tracks,
                   const CameraViewContainer& camera_views,
                   const KeysetContainer& keysets,
                   std::vector<Index>& keyset_id_offsets,
                   std::vector<size_t>& image_map,
                   std::vector<size_t>& point_map,
                   size_t& number_of_keys,
                   size_t& number_of_available_images,
                   size_t& number_of_available_points,
                   FeatureMapContainer& feature_maps) const
  {
    number_of_keys = 0;
    auto itr_keyset = keysets.begin();
    auto itr_keyset_end = keysets.end();
    keyset_id_offsets.resize(camera_views.size());
    for (size_t i = 0; itr_keyset != itr_keyset_end; ++itr_keyset, i++)
    {
      keyset_id_offsets[i] = number_of_keys;
      number_of_keys += itr_keyset->size();
    }

    image_map.resize(camera_views.size(), 0);
    auto itr_camera_view = camera_views.begin();
    auto itr_camera_view_end = camera_views.end();
    number_of_available_images = 0;
    for (size_t i = 0;
         itr_camera_view != itr_camera_view_end;
         ++itr_camera_view, ++i)
    {
      if (!itr_camera_view->empty())
      {
        image_map[i] = ++number_of_available_images;
      }
    }

    point_map.resize(tracks.size(), 0);
    auto itr_track = tracks.begin();
    auto itr_track_end = tracks.end();
    number_of_available_points = 0;
    for (size_t i = 0; itr_track != itr_track_end; ++itr_track, ++i)
    {
      if (!itr_track->empty())
      {
        point_map[i] = ++number_of_available_points;
      }
    }

    itr_track = tracks.begin();
    itr_track_end = tracks.end();
    feature_maps.resize(number_of_keys);
    for (size_t i = 0; itr_track != itr_track_end; ++itr_track, ++i)
    {
      auto itr_view = itr_track->begin();
      auto itr_view_end = itr_track->end();
      for (; itr_view != itr_view_end; ++itr_view)
      {
        Index feature_id = keyset_id_offsets[itr_view->first] +
                           itr_view->second;
        feature_maps[feature_id].first = Index(image_map[itr_view->first] -1);
        feature_maps[feature_id].second = Index(point_map[i] - 1);
      }
    }

    return 0;
  }

  Err GenerateVectorFunction(
    size_t number_of_available_images,
    size_t number_of_available_points,
    size_t number_of_keys,
    const std::vector<size_t>& image_map,
    const std::vector<size_t>& point_map,
    const TrackContainer& tracks,
    const CameraViewContainer& camera_views,
    const FeatureMapContainer& feature_maps,
    const ImageCameraMap& image_camera_map,
    const Point3DContainer& points,
    const ExtrinsicParamsContainer& extrinsic_params_set,
    std::vector<size_t>& constrained_image_ids,
    std::vector<size_t>& constrained_point_ids,
    VectorFunction& vector_function) const
  {
    vector_function.set_number_of_images(Index(number_of_available_images));
    vector_function.set_number_of_points(Index(number_of_available_points));
    vector_function.set_number_of_keys(Index(number_of_keys));
    vector_function.set_number_of_cameras(Index(intrinsic_params_set_.size()));
    vector_function.set_feature_maps(feature_maps);
    vector_function.set_image_camera_map(image_camera_map);
    //TODO:目前是计算所有的内参数，应根据mask选择性计算畸变和内参数
    vector_function.intrinsic_computations_mask().set();
    if (fix_mask_[FIX_POINTS])
    {
      //生成固定点
      BAPointContainer ba_points(number_of_available_points);
      auto itr_point = points.begin();
      auto itr_point_end = points.end();
      auto itr_track = tracks.begin();
      for (size_t i = 0; itr_point != itr_point_end; ++itr_point, ++itr_track)
      {
        if (!itr_track->empty())
        {
          ba_points[i] = *itr_point;
          i++;
        }
      }
      vector_function.set_fix_points(ba_points);
    }
    else
    {
      //生成点约束
      PointConstraintContainer point_constraints;
      for (size_t i = 0; i < number_of_planar_constrained_points_; i++)
      {
        while (1)
        {
          size_t point_id;
          IntegerUniformRandomGenerator<size_t>::Generate(
            0, tracks.size() - 1, point_id);
          auto itr_point = constrained_point_ids.begin();
          auto itr_point_end = constrained_point_ids.end();
          bool existed = false;
          for (; itr_point != itr_point_end; ++itr_point)
          {
            if (*itr_point == point_id)
            {
              existed = true;
              break;
            }
          }
          if (!existed && !tracks[point_id].empty())
          {
            PointConstraint point_constraint;
            point_constraint.point_id = point_map[point_id] - 1;
            point_constraint.mask.set(POINT_CONSTRAIN_X);
            point_constraint.mask.set(POINT_CONSTRAIN_Y);
            point_constraints.push_back(point_constraint);
            constrained_point_ids.push_back(point_id);
            break;
          }
        }
      }// for (size_t i = 0; i < number_of_planar_constrained_points_; i++)
      for (size_t i = 0; i < number_of_full_constrained_points_; i++)
      {
        while (1)
        {
          size_t point_id;
          IntegerUniformRandomGenerator<size_t>::Generate(
            0, tracks.size() - 1, point_id);
          auto itr_point = constrained_point_ids.begin();
          auto itr_point_end = constrained_point_ids.end();
          bool existed = false;
          for (; itr_point != itr_point_end; ++itr_point)
          {
            if (*itr_point == point_id)
            {
              existed = true;
              break;
            }
          }
          if (!existed && !tracks[point_id].empty())
          {
            PointConstraint point_constraint;
            point_constraint.point_id = point_map[point_id] - 1;
            point_constraint.mask.set(POINT_CONSTRAIN_X);
            point_constraint.mask.set(POINT_CONSTRAIN_Y);
            point_constraint.mask.set(POINT_CONSTRAIN_Z);
            point_constraints.push_back(point_constraint);
            constrained_point_ids.push_back(point_id);
            break;
          }
        }
      }// for (size_t i = 0; i < number_of_full_constrained_points_; i++)
      vector_function.point_constraints() = point_constraints;
    }

    if (fix_mask_[FIX_IMAGES])
    {
      //生成固定影像
      BAImageContainer ba_images(number_of_available_images);
      auto itr_extrin = extrinsic_params_set.begin();
      auto itr_extrin_end = extrinsic_params_set.end();
      auto itr_camera_view = camera_views.begin();
      for (size_t i = 0;
           itr_extrin != itr_extrin_end;
           ++itr_extrin, ++itr_camera_view)
      {
        if (!itr_camera_view->empty())
        {
          Vector3 t = -(itr_extrin->rotation() * itr_extrin->position());
          ba_images[i][0] = itr_extrin->rotation()[0];
          ba_images[i][1] = itr_extrin->rotation()[1];
          ba_images[i][2] = itr_extrin->rotation()[2];
          ba_images[i][3] = t[0];
          ba_images[i][4] = t[1];
          ba_images[i][5] = t[2];
          i++;
        }
      }
      vector_function.set_fix_images(ba_images);
    }
    else
    {
      //生成影像约束
      ImageConstraintContainer image_constraints;
      for (size_t i = 0; i < number_of_constrained_images_; i++)
      {
        while (1)
        {
          size_t image_id;
          IntegerUniformRandomGenerator<size_t>::Generate(
            0, camera_views.size() - 1, image_id);
          auto itr_image = constrained_image_ids.begin();
          auto itr_image_end = constrained_image_ids.end();
          bool existed = false;
          for (; itr_image != itr_image_end; ++itr_image)
          {
            if (*itr_image == image_id)
            {
              existed = true;
              break;
            }
          }
          if (!existed)
          {
            ImageConstraint image_constraint;
            image_constraint.image_id = image_map[image_id] - 1;
            image_constraint.mask.set();
            image_constraints.push_back(image_constraint);
            constrained_image_ids.push_back(image_id);
            break;
          }
        }
      }// for (size_t i = 0; i < number_of_constrained_images_; i++)
      vector_function.image_constraints() = image_constraints;
    }

    if (fix_mask_[FIX_CAMERAS])
    {
      //生成固定相机
      Index intrinsic_params_per_camera =
        vector_function.GetIntrinsicParamsSizePerCamera();
      BACameraContainer ba_cameras(intrinsic_params_set_.size(),
                                   BACamera(intrinsic_params_per_camera));
      for (size_t i = 0; i < intrinsic_params_set_.size(); i++)
      {
        //TODO:目前是计算所有的内参数，应根据mask选择性计算畸变和内参数
        ba_cameras[i][0] = intrinsic_params_set_[i].k1();
        ba_cameras[i][1] = intrinsic_params_set_[i].k2();
        ba_cameras[i][2] = intrinsic_params_set_[i].k3();
        ba_cameras[i][3] = intrinsic_params_set_[i].d1();
        ba_cameras[i][4] = intrinsic_params_set_[i].d2();
        ba_cameras[i][5] = intrinsic_params_set_[i].focal_length();
        ba_cameras[i][6] = intrinsic_params_set_[i].skew();
        ba_cameras[i][7] = intrinsic_params_set_[i].principal_point_x();
        ba_cameras[i][8] = intrinsic_params_set_[i].principal_point_y();
        ba_cameras[i][9] = intrinsic_params_set_[i].pixel_ratio();
      }
      vector_function.set_fix_cameras(ba_cameras);
    }
    else
    {
      //生成相机约束
      CameraConstraintContainer camera_constraints;
      for (size_t i = 0; i < number_of_constrained_cameras_; i++)
      {
        while (1)
        {
          size_t camera_id;
          IntegerUniformRandomGenerator<size_t>::Generate(
            0, intrinsic_params_set_.size() - 1, camera_id);
          auto itr_camera = camera_constraints.begin();
          auto itr_camera_end = camera_constraints.end();
          bool existed = false;
          for (; itr_camera != itr_camera_end; ++itr_camera)
          {
            if (itr_camera->camera_id == camera_id)
            {
              existed = true;
              break;
            }
          }
          if (!existed)
          {
            CameraConstraint camera_constraint;
            camera_constraint.camera_id = camera_id;
            camera_constraint.radial_mask.set();
            camera_constraint.decentering_mask.set();
            camera_constraint.intrinsic_mask.set();
            //camera_constraint.radial_mask.reset();
            //camera_constraint.decentering_mask.reset();
            //camera_constraint.intrinsic_mask.reset();
            //camera_constraint.intrinsic_mask.set(INTRINSIC_CONSTRAIN_SKEW);
            //camera_constraint.intrinsic_mask.set(
            //  INTRINSIC_CONSTRAIN_PIXEL_RATIO);
            camera_constraints.push_back(camera_constraint);
            break;
          }
        }
      }// for (size_t i = 0; i < number_of_constrained_cameras_; i++)
      vector_function.camera_constraints() = camera_constraints;
    }

    return 0;
  }

  Err GenerateXVector(const Point3DContainer& points,
                      const TrackContainer& tracks,
                      const ExtrinsicParamsContainer& extrinsic_params_set,
                      const CameraViewContainer& camera_views,
                      const VectorFunction& vector_function,
                      XVector& x) const
  {
    Index x_size = vector_function.GetXSize();
    x.resize(x_size);

    if (!fix_mask_[FIX_POINTS])
    {
      auto itr_point = points.begin();
      auto itr_point_end = points.end();
      auto itr_track = tracks.begin();
      auto itr_track_end = tracks.end();
      for (Index i = 0; itr_point != itr_point_end; ++itr_point, ++itr_track)
      {
        if (!itr_track->empty())
        {
          x.segment(i * VectorFunction::params_per_point_,
                    VectorFunction::params_per_point_) = *itr_point;
          ++i;
        }
      }
    }

    Index x_offset =
      fix_mask_[FIX_POINTS] ? 0 : vector_function.GetPointParamsSize();
    if (!fix_mask_[FIX_IMAGES])
    {
      auto itr_extrin = extrinsic_params_set.begin();
      auto itr_extrin_end = extrinsic_params_set.end();
      auto itr_camera_view = camera_views.begin();
      auto itr_camera_view_end = camera_views.end();
      for (Index i = 0;
           itr_extrin != itr_extrin_end;
           ++itr_extrin, ++itr_camera_view)
      {
        if (!itr_camera_view->empty())
        {
          Vector3 t = -(itr_extrin->rotation() * itr_extrin->position());
          x[x_offset + i * VectorFunction::extrinsic_params_per_image_ + 0] =
            itr_extrin->rotation()[0];
          x[x_offset + i * VectorFunction::extrinsic_params_per_image_ + 1] =
            itr_extrin->rotation()[1];
          x[x_offset + i * VectorFunction::extrinsic_params_per_image_ + 2] =
            itr_extrin->rotation()[2];
          x[x_offset + i * VectorFunction::extrinsic_params_per_image_ + 3] =
            t[0];
          x[x_offset + i * VectorFunction::extrinsic_params_per_image_ + 4] =
            t[1];
          x[x_offset + i * VectorFunction::extrinsic_params_per_image_ + 5] =
            t[2];

          i++;
        }
      }
    }

    x_offset +=
      fix_mask_[FIX_IMAGES] ? 0 : vector_function.GetExtrinsicParamsSize();
    if (!fix_mask_[FIX_CAMERAS])
    {
      Index intrinsic_params_size_per_camera =
        vector_function.GetIntrinsicParamsSizePerCamera();
      for (Index i = 0; i < Index(intrinsic_params_set_.size()); i++)
      {
        const IntrinsicParams& intrinsic_params = intrinsic_params_set_[i];
        //TODO:目前是计算所有的内参数，应根据mask选择性计算畸变和内参数
        x[x_offset + i * intrinsic_params_size_per_camera + 0] =
          intrinsic_params.k1();
        x[x_offset + i * intrinsic_params_size_per_camera + 1] =
          intrinsic_params.k2();
        x[x_offset + i * intrinsic_params_size_per_camera + 2] =
          intrinsic_params.k3();
        x[x_offset + i * intrinsic_params_size_per_camera + 3] =
          intrinsic_params.d1();
        x[x_offset + i * intrinsic_params_size_per_camera + 4] =
          intrinsic_params.d2();
        x[x_offset + i * intrinsic_params_size_per_camera + 5] =
          intrinsic_params.focal_length();
        x[x_offset + i * intrinsic_params_size_per_camera + 6] =
          intrinsic_params.skew();
        x[x_offset + i * intrinsic_params_size_per_camera + 7] =
          intrinsic_params.principal_point_x();
        x[x_offset + i * intrinsic_params_size_per_camera + 8] =
          intrinsic_params.principal_point_y();
        x[x_offset + i * intrinsic_params_size_per_camera + 9] =
          intrinsic_params.pixel_ratio();
      }
    }

    return 0;
  }

  Err GenerateYVector(const KeysetContainer& keysets,
                      const ExtrinsicParamsContainer& extrinsic_params_set,
                      const Point3DContainer& points,
                      const TrackContainer& tracks,
                      const VectorFunction& vector_function,
                      const std::vector<size_t>& constrained_image_ids,
                      const std::vector<size_t>& constrained_point_ids,
                      size_t number_of_available_points,
                      YVector& y) const
  {
    Index y_size = vector_function.GetYSize();
    y.resize(y_size);

    auto itr_keyset = keysets.begin();
    auto itr_keyset_end = keysets.end();
    for (Index i = 0; itr_keyset != itr_keyset_end; ++itr_keyset)
    {
      for (size_t key_id = 0; key_id < itr_keyset->size(); key_id++)
      {
        //TODO:加入选择性计算畸变和内参数的mask，并根据该mask转换合适的key
        y.segment(i * VectorFunction::params_per_key_,
                  VectorFunction::params_per_key_) =
          (*itr_keyset)[key_id];
        i++;
      }
    }

    Index y_offset = vector_function.GetYKeysSize();

    if (!fix_mask_[FIX_POINTS])
    {
      auto itr_point_id = constrained_point_ids.begin();
      auto itr_point_id_end = constrained_point_ids.end();
      for (size_t i = 0; itr_point_id != itr_point_id_end; ++itr_point_id, ++i)
      {
        const PointConstraint& point_constraint =
          vector_function.point_constraints()[i];
        if (point_constraint.mask[POINT_CONSTRAIN_X])
        {
          y[y_offset] = points[*itr_point_id][0];
          y_offset++;
        }
        if (point_constraint.mask[POINT_CONSTRAIN_Y])
        {
          y[y_offset] = points[*itr_point_id][1];
          y_offset++;
        }
        if (point_constraint.mask[POINT_CONSTRAIN_Z])
        {
          y[y_offset] = points[*itr_point_id][2];
          y_offset++;
        }
      }
    }

    if (!fix_mask_[FIX_IMAGES])
    {
      auto itr_image_id = constrained_image_ids.begin();
      auto itr_image_id_end = constrained_image_ids.end();
      for (size_t i = 0; itr_image_id != itr_image_id_end; ++itr_image_id, ++i)
      {
        const ImageConstraint& image_constraint =
          vector_function.image_constraints()[i];
        const ExtrinsicParams& extrinsic_params =
          extrinsic_params_set[*itr_image_id];
        Vector3 t =
          -(extrinsic_params.rotation() * extrinsic_params.position());
        if (image_constraint.mask[IMAGE_CONSTRAIN_ROTATION])
        {
          y[y_offset + 0] = extrinsic_params.rotation()[0];
          y[y_offset + 1] = extrinsic_params.rotation()[1];
          y[y_offset + 2] = extrinsic_params.rotation()[2];
          y_offset += 3;
        }
        if (image_constraint.mask[IMAGE_CONSTRAIN_POSITION_X])
        {
          y[y_offset] = t[0];
          y_offset++;
        }
        if (image_constraint.mask[IMAGE_CONSTRAIN_POSITION_Y])
        {
          y[y_offset] = t[1];
          y_offset++;
        }
        if (image_constraint.mask[IMAGE_CONSTRAIN_POSITION_Z])
        {
          y[y_offset] = t[2];
          y_offset++;
        }
      }
    }

    if (!fix_mask_[FIX_CAMERAS])
    {
      for (size_t i = 0; i < number_of_constrained_cameras_; i++)
      {
        const CameraConstraint& camera_constraint =
          vector_function.camera_constraints()[i];
        size_t camera_id = camera_constraint.camera_id;
        if (camera_constraint.radial_mask[RADIAL_CONSTRAIN_K1])
        {
          y[y_offset] = intrinsic_params_set_[camera_id].k1();
          y_offset++;
        }
        if (camera_constraint.radial_mask[RADIAL_CONSTRAIN_K2])
        {
          y[y_offset] = intrinsic_params_set_[camera_id].k2();
          y_offset++;
        }
        if (camera_constraint.radial_mask[RADIAL_CONSTRAIN_K3])
        {
          y[y_offset] = intrinsic_params_set_[camera_id].k3();
          y_offset++;
        }
        if (camera_constraint.decentering_mask[DECENTERING_CONSTRAIN_D1])
        {
          y[y_offset] = intrinsic_params_set_[camera_id].d1();
          y_offset++;
        }
        if (camera_constraint.decentering_mask[DECENTERING_CONSTRAIN_D2])
        {
          y[y_offset] = intrinsic_params_set_[camera_id].d2();
          y_offset++;
        }
        if (camera_constraint.intrinsic_mask[INTRINSIC_CONSTRAIN_FOCAL_LENGTH])
        {
          y[y_offset] = intrinsic_params_set_[camera_id].focal_length();
          y_offset++;
        }
        if (camera_constraint.intrinsic_mask[INTRINSIC_CONSTRAIN_SKEW])
        {
          y[y_offset] = intrinsic_params_set_[camera_id].skew();
          y_offset++;
        }
        if (camera_constraint.intrinsic_mask[INTRINSIC_CONSTRAIN_PRINCIPAL_X])
        {
          y[y_offset] = intrinsic_params_set_[camera_id].principal_point_x();
          y_offset++;
        }
        if (camera_constraint.intrinsic_mask[INTRINSIC_CONSTRAIN_PRINCIPAL_Y])
        {
          y[y_offset] = intrinsic_params_set_[camera_id].principal_point_y();
          y_offset++;
        }
        if (camera_constraint.intrinsic_mask[INTRINSIC_CONSTRAIN_PIXEL_RATIO])
        {
          y[y_offset] = intrinsic_params_set_[camera_id].pixel_ratio();
          y_offset++;
        }
      }// for (size_t i = 0; i < number_of_constrained_cameras_; i++)
    }

    return 0;
  }

private:
  MultipleFlightGenerator multiple_flight_generator_;
  KeysetGenerator keyset_generator_;
  size_t number_of_points_;
  size_t number_of_planar_constrained_points_;
  size_t number_of_full_constrained_points_;
  size_t number_of_constrained_images_;
  size_t number_of_constrained_cameras_;
  IntrinsicParamsContainer intrinsic_params_set_;
  FixMask fix_mask_;
};

}
}
}

#endif
