#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "gtest/gtest.h"

#include "hs_sfm/synthetic/multiple_flight_generator.hpp"
#include "hs_sfm/synthetic/multiple_camera_keyset_generator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestMultipleCameraGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

public:
  typedef hs::sfm::synthetic::MultipleFlightGenerator<Scalar, ImageDimension>
          MultipleFlightGenerator;
  typedef typename MultipleFlightGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename MultipleFlightGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename MultipleFlightGenerator::Image Image;
  typedef typename MultipleFlightGenerator::ImageContainer ImageContainer;
  typedef typename MultipleFlightGenerator::Point3D Point3D;
  typedef typename MultipleFlightGenerator::Point3DContainer Point3DContainer;
  typedef typename MultipleFlightGenerator::FlightGenerator FlightGenerator;
  typedef typename MultipleFlightGenerator::FlightGeneratorContainer
                   FlightGeneratorContainer;

  typedef hs::sfm::synthetic::MultipleCameraKeysetGenerator<Scalar,
                                                            ImageDimension>
          KeysetGenerator;
  typedef typename KeysetGenerator::IntrinsicParams IntrinsicParams;
  typedef typename KeysetGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

private:
  typedef typename Keyset::Key Key;
  typedef hs::sfm::ProjectiveFunctions<Scalar> ProjectiveFunctions;

public:
  TestMultipleCameraGenerator(
    Scalar flight_longitudinal_overlap_ratio,
    Scalar flight_lateral_overlap_ratio,
    Scalar north_west_angle,
    Scalar north_west_angle_stddev,
    Scalar offset_stddev,
    const FlightGeneratorContainer& flight_generators,
    const IntrinsicParamsContainer& intrinsic_params_set,
    size_t number_of_points,
    const std::string& test_name)
    : multiple_flight_generator_(flight_longitudinal_overlap_ratio,
                                 flight_lateral_overlap_ratio,
                                 north_west_angle,
                                 north_west_angle_stddev,
                                 offset_stddev,
                                 flight_generators),
      keyset_generator_(),
      intrinsic_params_set_(intrinsic_params_set),
      number_of_points_(number_of_points),
      test_name_(test_name) {}

  Err Test() const
  {
    Point3DContainer points;
    ExtrinsicParamsContainer extrinsic_params_set;
    ImageContainer images;
    std::vector<size_t> image_intrinsic_map;
    KeysetContainer keysets;
    hs::sfm::TrackContainer tracks;
    hs::sfm::CameraViewContainer camera_views;
    if (TestGenerator(points,
                      extrinsic_params_set,
                      images,
                      image_intrinsic_map,
                      keysets,
                      tracks,
                      camera_views) != 0)
    {
      std::cout<<"Testing generator failed!\n";
      return -1;
    }

    if (TestExtrinsicParamsSets(extrinsic_params_set, images) != 0)
    {
      std::cout<<"Testing extrinsic params sets failed!\n";
      return -1;
    }

    if (TestPoints(points) != 0)
    {
      std::cout<<"Testing points failed!\n";
      return -1;
    }

    if (TestReprojectiveError(extrinsic_params_set,
                              images,
                              image_intrinsic_map,
                              points,
                              keysets,
                              tracks,
                              camera_views) != 0)
    {
      std::cout<<"Testing reprojective error failed!\n";
      return -1;
    }

    return 0;
  }

private:
  Err TestGenerator(Point3DContainer& points,
                    ExtrinsicParamsContainer& extrinsic_params_set,
                    ImageContainer& images,
                    std::vector<size_t>& image_intrinsic_map,
                    KeysetContainer& keysets,
                    hs::sfm::TrackContainer& tracks,
                    hs::sfm::CameraViewContainer& camera_views) const
  {
    if (multiple_flight_generator_.GeneratePoints(number_of_points_,
                                                  points) != 0) return -1;

    size_t number_of_flights = multiple_flight_generator_.GetNumberOfFlights();
    extrinsic_params_set.clear();
    images.clear();

    for (size_t i = 0; i < number_of_flights; i++)
    {
      ExtrinsicParamsContainer extrinsic_params_set_flight;
      ImageContainer images_flight;
      if (multiple_flight_generator_.GenerateExtrinsicParamsContainer(
            i, extrinsic_params_set_flight, images_flight) != 0) return -1;

      for (size_t j = 0; j < images_flight.size(); j++)
      {
        extrinsic_params_set.push_back(extrinsic_params_set_flight[j]);
        images.push_back(images_flight[j]);
        image_intrinsic_map.push_back(i);
      }
    }

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

  Err TestExtrinsicParamsSets(
    const ExtrinsicParamsContainer& extrinsic_params_set,
    const ImageContainer& images) const
  {
    size_t number_of_images = images.size();
    if (number_of_images != extrinsic_params_set.size())
    {
      return -1;
    }

    std::string path = test_name_ + "_extrinisc_params.txt";
    std::ofstream file(path.c_str(), std::ios::out);
    if (!file.is_open())
    {
      return -1;
    }
    for (size_t i = 0; i < number_of_images; i++)
    {
      file<<i<<" "
          <<extrinsic_params_set[i].position()[0]<<" "
          <<extrinsic_params_set[i].position()[1]<<" "
          <<extrinsic_params_set[i].position()[2]<<"\n";
    }

    return 0;
  }

  Err TestPoints(const Point3DContainer& points) const
  {
    if (points.size() != number_of_points_)
    {
      return -1;
    }

    std::string path = test_name_ + "_points.txt";
    std::ofstream file(path.c_str(), std::ios::out);
    if (!file.is_open())
    {
      return -1;
    }
    for (size_t i = 0; i < number_of_points_; i++)
    {
      file<<i<<" "
          <<points[i][0]<<" "
          <<points[i][1]<<" "
          <<points[i][2]<<"\n";
    }

    return 0;
  }

  Err TestReprojectiveError(
    const ExtrinsicParamsContainer& extrinsic_params_set,
    const ImageContainer& images,
    const std::vector<size_t>& image_intrinsic_map,
    const Point3DContainer& points,
    const KeysetContainer& keysets,
    const hs::sfm::TrackContainer& tracks,
    const hs::sfm::CameraViewContainer& camera_views) const
  {
    size_t number_of_tracks = tracks.size();
    Scalar threshold = Scalar(1e-8);

    size_t number_of_track_reprojections = 0;
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();
      const Point3D& point = points[i];
      for (size_t j = 0; j < number_of_views; j++)
      {
        size_t image_id = tracks[i][j].first;
        size_t key_id = tracks[i][j].second;
        size_t camera_id = image_intrinsic_map[image_id];
        const IntrinsicParams& intrinsic_params =
          intrinsic_params_set_[camera_id];
        const ExtrinsicParams& extrinsic_params =
          extrinsic_params_set[image_id];
        Key key = ProjectiveFunctions::WorldPointProjectToImageKey(
                    intrinsic_params,
                    extrinsic_params,
                    point);
        if (!key.isApprox(keysets[image_id][key_id], threshold))
        {
          return -1;
        }

        number_of_track_reprojections++;
      }// for (size_t j = 0; j < number_of_views; j++)
    }// for (size_t i = 0; i < number_of_tracks; i++)

    size_t number_of_images = images.size();
    size_t number_of_camera_views_reprojections = 0;
    for (size_t i = 0; i < number_of_images; i++)
    {
      size_t camera_id = image_intrinsic_map[i];
      const IntrinsicParams& intrinsic_params =
        intrinsic_params_set_[camera_id];
      const ExtrinsicParams& extrinsic_params =
        extrinsic_params_set[i];
      size_t number_of_views = camera_views[i].size();
      for (size_t j = 0; j < number_of_views; j++)
      {
        size_t point_id = camera_views[i][j].first;
        size_t key_id = camera_views[i][j].second;
        const Point3D& point = points[point_id];
        Key key = ProjectiveFunctions::WorldPointProjectToImageKey(
                    intrinsic_params,
                    extrinsic_params,
                    point);
        if (!key.isApprox(keysets[i][key_id], threshold))
        {
          return -1;
        }

        number_of_camera_views_reprojections++;
      }// for (size_t j = 0; j < number_of_views; j++)
    }// for (size_t i = 0; i < number_of_images; i++)

    if (number_of_track_reprojections != number_of_camera_views_reprojections)
    {
      return -1;
    }

    return 0;
  }

private:
  MultipleFlightGenerator multiple_flight_generator_;
  KeysetGenerator keyset_generator_;
  IntrinsicParamsContainer intrinsic_params_set_;
  size_t number_of_points_;
  std::string test_name_;
};

TEST(TestMultipleCameraGenerator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestMultipleCameraGenerator<Scalar, ImageDimension> Test;
  typedef Test::FlightGenerator FlightGenerator;
  typedef Test::FlightGeneratorContainer FlightGeneratorContainer;
  typedef Test::IntrinsicParams IntrinsicParams;
  typedef Test::IntrinsicParamsContainer IntrinsicParamsContainer;

  FlightGeneratorContainer flight_generators;
  IntrinsicParamsContainer intrinsic_params_set;

  Scalar focal_length_in_metre_0 = 0.018858358970276164;
  size_t number_of_strips_0 = 15;
  size_t number_of_cameras_in_strip_0 = 20;
  Scalar ground_resolution_0 = 0.1;
  ImageDimension image_width_0 = 6000;
  ImageDimension image_height_0 = 4000;
  Scalar pixel_size_0 = 0.0000039;
  size_t number_of_points_0 = 2000;
  Scalar lateral_overlap_ratio_0 = 0.6;
  Scalar longitudinal_overlap_ratio_0 = 0.8;
  Scalar scene_max_height_0 = 50;
  Scalar camera_height_stddev_0 = 2;
  Scalar camera_planar_stddev_0 = 2;
  Scalar camera_rotation_stddev_0 = 10;
  FlightGenerator flight_generator_0(
    focal_length_in_metre_0,
    number_of_strips_0,
    number_of_cameras_in_strip_0,
    ground_resolution_0,
    image_width_0,
    image_height_0,
    pixel_size_0,
    number_of_points_0,
    lateral_overlap_ratio_0,
    longitudinal_overlap_ratio_0,
    scene_max_height_0,
    camera_height_stddev_0,
    camera_planar_stddev_0,
    camera_rotation_stddev_0);
  flight_generators.push_back(flight_generator_0);
  IntrinsicParams intrinsic_params_0(focal_length_in_metre_0 / pixel_size_0,
                                     0,
                                     -42.4095312016,
                                     -31.699212823,
                                     1,
                                     -0.0050490462006048831,
                                     0.031293804298609881,
                                     -0.030794820960442223,
                                     -0.00055376548320189127,
                                     -0.00049877717768381476);
  intrinsic_params_set.push_back(intrinsic_params_0);

  Scalar focal_length_in_metre_1 = 0.02995452167701055;
  size_t number_of_strips_1 = 20;
  size_t number_of_cameras_in_strip_1 = 15;
  Scalar ground_resolution_1 = 0.1;
  ImageDimension image_width_1 = 6000;
  ImageDimension image_height_1 = 4000;
  Scalar pixel_size_1 = 0.0000039;
  size_t number_of_points_1 = 2000;
  Scalar lateral_overlap_ratio_1 = 0.6;
  Scalar longitudinal_overlap_ratio_1 = 0.8;
  Scalar scene_max_height_1 = 50;
  Scalar camera_height_stddev_1 = 2;
  Scalar camera_planar_stddev_1 = 2;
  Scalar camera_rotation_stddev_1 = 10;
  FlightGenerator flight_generator_1(
    focal_length_in_metre_1,
    number_of_strips_1,
    number_of_cameras_in_strip_1,
    ground_resolution_1,
    image_width_1,
    image_height_1,
    pixel_size_1,
    number_of_points_1,
    lateral_overlap_ratio_1,
    longitudinal_overlap_ratio_1,
    scene_max_height_1,
    camera_height_stddev_1,
    camera_planar_stddev_1,
    camera_rotation_stddev_1);
  flight_generators.push_back(flight_generator_1);
  IntrinsicParams intrinsic_params_1(focal_length_in_metre_1 / pixel_size_1,
                                     0,
                                     -21.669436058,
                                     -44.8644764322,
                                     1,
                                     -0.02529179096221609,
                                     0.23762413973445157,
                                     -0.64208397668697237,
                                     -0.0020605099808780948,
                                     -0.00028706423764766859);
  intrinsic_params_set.push_back(intrinsic_params_1);

  Scalar focal_length_in_metre_2 = 0.019056097774998712;
  size_t number_of_strips_2 = 10;
  size_t number_of_cameras_in_strip_2 = 30;
  Scalar ground_resolution_2 = 0.1;
  ImageDimension image_width_2 = 6000;
  ImageDimension image_height_2 = 4000;
  Scalar pixel_size_2 = 0.0000039;
  size_t number_of_points_2 = 2000;
  Scalar lateral_overlap_ratio_2 = 0.6;
  Scalar longitudinal_overlap_ratio_2 = 0.8;
  Scalar scene_max_height_2 = 50;
  Scalar camera_height_stddev_2 = 2;
  Scalar camera_planar_stddev_2 = 2;
  Scalar camera_rotation_stddev_2 = 10;
  FlightGenerator flight_generator_2(
    focal_length_in_metre_2,
    number_of_strips_2,
    number_of_cameras_in_strip_2,
    ground_resolution_2,
    image_width_2,
    image_height_2,
    pixel_size_2,
    number_of_points_2,
    lateral_overlap_ratio_2,
    longitudinal_overlap_ratio_2,
    scene_max_height_2,
    camera_height_stddev_2,
    camera_planar_stddev_2,
    camera_rotation_stddev_2);
  flight_generators.push_back(flight_generator_2);
  IntrinsicParams intrinsic_params_2(focal_length_in_metre_2 / pixel_size_2,
                                     0,
                                     -35.2052431556,
                                     -16.4262220759,
                                     1,
                                     -0.10316088386868619,
                                     0.13520490482776426,
                                     -0.05489235547426094,
                                     4.1434720317373253e-006,
                                     -0.00025018439997095336);
  intrinsic_params_set.push_back(intrinsic_params_2);

  Test test(0.8, 0.2, 60, 10, 15, flight_generators,
            intrinsic_params_set,
            5000, "simple_test");
  ASSERT_EQ(0, test.Test());
}

}
