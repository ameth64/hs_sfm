#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_GENERAL_SYNTHETIC_DATA_GENERATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_GENERAL_SYNTHETIC_DATA_GENERATOR_HPP_

#include "hs_sfm/synthetic/multiple_flight_generator.hpp"
#include "hs_sfm/synthetic/multiple_camera_keyset_generator.hpp"
#include "hs_sfm/bundle_adjustment/general_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar, typename _ImageDimension>
class GeneralSyntheticDataGenerator
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

  typedef GeneralVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::FeatureMap FeatureMap;
  typedef typename VectorFunction::FeatureMapContainer FeatureMapContainer;
  typedef typename VectorFunction::ImageCameraMap ImageCameraMap;
  typedef typename VectorFunction::IntrinsicConstraintsMask
                   IntrinsicConstraintsMask;
  typedef typename VectorFunction::StructureConstraintsMask
                   StructureConstraintsMask;

private:
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;

public:
  GeneralSyntheticDataGenerator(
    Scalar flight_longitudinal_overlap_ratio,
    Scalar flight_lateral_overlap_ratio,
    Scalar north_west_angle,
    Scalar north_west_angle_stddev,
    Scalar offset_stddev,
    const FlightGeneratorContainer& flight_generators,
    size_t number_of_points,
    size_t number_of_constrained_points,
    const IntrinsicParamsContainer& intrinsic_params_set,
    const IntrinsicConstraintsMask& intrinsic_constraints_mask,
    const StructureConstraintsMask& structure_constraints_mask)
    : multiple_flight_generator_(flight_longitudinal_overlap_ratio,
                                 flight_lateral_overlap_ratio,
                                 north_west_angle,
                                 north_west_angle_stddev,
                                 offset_stddev,
                                 flight_generators),
      keyset_generator_(),
      number_of_points_(number_of_points),
      number_of_constrained_points_(number_of_constrained_points),
      intrinsic_params_set_(intrinsic_params_set),
      intrinsic_constraints_mask_(intrinsic_constraints_mask),
      structure_constraints_mask_(structure_constraints_mask) {}

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

    if (GenerateVectorFunction(number_of_available_images,
                               number_of_available_points,
                               number_of_keys,
                               feature_maps,
                               image_camera_map,
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
                        points,
                        tracks,
                        vector_function,
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

  Err GenerateVectorFunction(size_t number_of_available_images,
                             size_t number_of_available_points,
                             size_t number_of_keys,
                             const FeatureMapContainer& feature_maps,
                             const ImageCameraMap& image_camera_map,
                             VectorFunction& vector_function) const
  {
    vector_function.set_number_of_images(Index(number_of_available_images));
    vector_function.set_number_of_points(Index(number_of_available_points));
    vector_function.set_number_of_keys(Index(number_of_keys));
    vector_function.set_number_of_constrained_points(
      Index(number_of_constrained_points_));
    vector_function.set_number_of_cameras(Index(intrinsic_params_set_.size()));
    vector_function.set_feature_maps(feature_maps);
    vector_function.set_image_camera_map(image_camera_map);
    //TODO:目前是计算所有的内参数，应根据mask选择性计算畸变和内参数
    vector_function.intrinsic_computations_mask().set();
    vector_function.intrinsic_constraints_mask() = intrinsic_constraints_mask_;
    vector_function.structure_constraints_mask() = structure_constraints_mask_;
    return -1;
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

    Index x_offset = vector_function.GetPointParamsSize();
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

    x_offset += vector_function.GetExtrinsicParamsSize();
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

    return 0;
  }

  Err GenerateYVector(const KeysetContainer& keysets,
                      const Point3DContainer& points,
                      const TrackContainer& tracks,
                      const VectorFunction& vector_function,
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
    if (intrinsic_constraints_mask_.any())
    {
      for (Index i = 0; i < Index(intrinsic_params_set_.size()); i++)
      {
        const IntrinsicParams& intrinsic_params = intrinsic_params_set_[i];
        Index y_intrinsic_offset =
          y_offset + i * Index(intrinsic_constraints_mask_.count());
        if (intrinsic_constraints_mask_[
              VectorFunction::CONSTRAIN_RADIAL_K1])
        {
          y[y_intrinsic_offset] = intrinsic_params.k1();
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[
              VectorFunction::CONSTRAIN_RADIAL_K2])
        {
          y[y_intrinsic_offset] = intrinsic_params.k2();
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[
              VectorFunction::CONSTRAIN_RADIAL_K3])
        {
          y[y_intrinsic_offset] = intrinsic_params.k3();
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[
              VectorFunction::CONSTRAIN_DECENTERING_D1])
        {
          y[y_intrinsic_offset] = intrinsic_params.d1();
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[
              VectorFunction::CONSTRAIN_DECENTERING_D2])
        {
          y[y_intrinsic_offset] = intrinsic_params.d2();
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[
              VectorFunction::CONSTRAIN_FOCAL_LENGTH])
        {
          y[y_intrinsic_offset] = intrinsic_params.focal_length();
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[
              VectorFunction::CONSTRAIN_SKEW])
        {
          y[y_intrinsic_offset] = intrinsic_params.skew();
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[
              VectorFunction::CONSTRAIN_PRINCIPAL_X])
        {
          y[y_intrinsic_offset] = intrinsic_params.principal_point_x();
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[
              VectorFunction::CONSTRAIN_PRINCIPAL_Y])
        {
          y[y_intrinsic_offset] = intrinsic_params.principal_point_y();
          y_intrinsic_offset++;
        }
        if (intrinsic_constraints_mask_[
              VectorFunction::CONSTRAIN_PIXEL_RATIO])
        {
          y[y_intrinsic_offset] = intrinsic_params.pixel_ratio();
          y_intrinsic_offset++;
        }

      }// for (Index i = 0; i < Index(intrinsic_params_set_.size()); i++)
    }// if (intrinsic_constraints_mask_.any())

    y_offset += vector_function.GetYConstrainedIntrinsicParamsSize();
    if (structure_constraints_mask_.any())
    {
      if (structure_constraints_mask_[VectorFunction::CONSTRAIN_POINTS])
      {
        size_t constrained_points_begin = number_of_available_points -
                                          number_of_constrained_points_;
        auto itr_point = points.begin();
        auto itr_point_end = points.end();
        auto itr_track = tracks.begin();
        auto itr_track_end = tracks.end();
        for (size_t i = 0; itr_point != itr_point_end; ++itr_point, ++itr_track)
        {
          if (!itr_track->empty())
          {
            if (i >= constrained_points_begin)
            {
              Index constrained_point_id = Index(i - constrained_points_begin);
              y.segment(y_offset +
                        constrained_point_id *
                        VectorFunction::params_per_point_,
                        VectorFunction::params_per_point_) = *itr_point;
            }
            i++;
          }
        }
      }// if (structure_constraints_mask_[CONSTRAIN_POINTS])
    }// if (structure_constraints_mask_.any())

    return 0;
  }

private:
  MultipleFlightGenerator multiple_flight_generator_;
  KeysetGenerator keyset_generator_;
  size_t number_of_points_;
  size_t number_of_constrained_points_;
  IntrinsicParamsContainer intrinsic_params_set_;
  IntrinsicConstraintsMask intrinsic_constraints_mask_;
  StructureConstraintsMask structure_constraints_mask_;

};

}
}
}

#endif
