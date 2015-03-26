#include <gtest/gtest.h>

#include "data_tester.hpp"
#include "synthetic_data_generator.hpp"
#include "real_data_generator.hpp"

#include "hs_test_utility/test_env/data_path.hpp"

#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"
#include "hs_sfm/incremental/incremental.hpp"
#include "hs_sfm/incremental/gcp_similar_transform_estimator.hpp"

namespace
{

template <typename _Scalar>
class TestIncremental
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

private:
  typedef hs::sfm::incremental::IncrementalSFM<Scalar> IncrementalSFM;
  typedef hs::sfm::ObjectIndexMap ObjectIndexMap;
  typedef hs::sfm::ViewInfoIndexer ViewInfoIndexer;
  typedef hs::sfm::incremental::DataTester<Scalar> Tester;
  typedef hs::sfm::incremental::GCPSimilarTransformEstimator<Scalar>
          SimilarTransformEstimator;
  typedef typename SimilarTransformEstimator::Rotation Rotation;
  typedef typename SimilarTransformEstimator::Translate Translate;

public:
  typedef typename IncrementalSFM::IntrinsicParams IntrinsicParams;
  typedef typename IncrementalSFM::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename IncrementalSFM::ExtrinsicParams ExtrinsicParams;
  typedef typename IncrementalSFM::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename IncrementalSFM::KeysetContainer KeysetContainer;
  typedef typename IncrementalSFM::Point Point;
  typedef typename IncrementalSFM::PointContainer PointContainer;
  typedef hs::sfm::MatchContainer MatchContainer;
  typedef hs::sfm::TrackContainer TrackContainer;

public:
  TestIncremental(
    Scalar ground_resolution,
    Scalar key_stddev,
    const std::string& test_name)
    : key_stddev_(key_stddev),
      test_name_(test_name),
      ground_resolution_(ground_resolution) {}

  Err operator() (
    const IntrinsicParamsContainer& intrinsic_params_set_initial,
    const std::vector<size_t>& image_intrinsic_map_input,
    const MatchContainer& matches,
    const KeysetContainer& keysets,
    const PointContainer& gcps,
    const TrackContainer& tracks_gcp,
    const KeysetContainer& keysets_gcp,
    const IntrinsicParamsContainer& intrinsic_params_set_true =
      IntrinsicParamsContainer(),
    const ExtrinsicParamsContainer& extrinsic_params_set_absolute_true =
      ExtrinsicParamsContainer(),
    const PointContainer& points_absolute_true =
      PointContainer()) const
  {
    if (keysets.size() != image_intrinsic_map_input.size()) return -1;
    ExtrinsicParamsContainer extrinsic_params_set_relative_estimate;
    ObjectIndexMap image_intrinsic_map(image_intrinsic_map_input.size());
    for (size_t i = 0; i < image_intrinsic_map_input.size(); i++)
    {
      image_intrinsic_map.SetObjectId(i, image_intrinsic_map_input[i]);
    }
    ObjectIndexMap image_extrinsic_map;
    PointContainer points_relative_estimate;
    TrackContainer tracks;
    ObjectIndexMap track_point_map;
    ViewInfoIndexer view_info_indexer;
    IncrementalSFM incremental_sfm(100, 8, 2, 7);
    IntrinsicParamsContainer intrinsic_params_set_estimate =
      intrinsic_params_set_initial;
    if (incremental_sfm(image_intrinsic_map,
                        matches,
                        keysets,
                        intrinsic_params_set_estimate,
                        extrinsic_params_set_relative_estimate,
                        image_extrinsic_map,
                        points_relative_estimate,
                        tracks,
                        track_point_map,
                        view_info_indexer) != 0)
    {
      return -1;
    }

    for (size_t i = 0; i < intrinsic_params_set_estimate.size(); i++)
    {
      std::cout<<"intrinsic_params_set_estimate["<<i<<"]:\n";
      std::cout<<intrinsic_params_set_estimate[i].focal_length()<<"\n";
      std::cout<<intrinsic_params_set_estimate[i].skew()<<"\n";
      std::cout<<intrinsic_params_set_estimate[i].principal_point_x()<<"\n";
      std::cout<<intrinsic_params_set_estimate[i].principal_point_y()<<"\n";
      std::cout<<intrinsic_params_set_estimate[i].pixel_ratio()<<"\n";
      std::cout<<intrinsic_params_set_estimate[i].k1()<<"\n";
      std::cout<<intrinsic_params_set_estimate[i].k2()<<"\n";
      std::cout<<intrinsic_params_set_estimate[i].k3()<<"\n";
      std::cout<<intrinsic_params_set_estimate[i].d1()<<"\n";
      std::cout<<intrinsic_params_set_estimate[i].d2()<<"\n";
    }

    Tester tester;
    if (tester.TestReprojectiveError(
          keysets,
          intrinsic_params_set_estimate,
          tracks,
          image_intrinsic_map,
          image_extrinsic_map,
          track_point_map,
          view_info_indexer,
          extrinsic_params_set_relative_estimate,
          points_relative_estimate,
          key_stddev_) != 0) return -1;

    SimilarTransformEstimator similar_transform_estimator;
    ObjectIndexMap track_point_map_gcp;
    ViewInfoIndexer view_info_indexer_gcp;
    PointContainer gcps_relative;
    Rotation rotation_similar;
    Translate translate_similar;
    Scalar scale_similar;
    if (similar_transform_estimator(keysets_gcp,
                                    intrinsic_params_set_estimate,
                                    extrinsic_params_set_relative_estimate,
                                    tracks_gcp,
                                    gcps,
                                    image_intrinsic_map,
                                    image_extrinsic_map,
                                    4,
                                    key_stddev_ * Scalar(4),
                                    rotation_similar,
                                    translate_similar,
                                    scale_similar,
                                    track_point_map_gcp,
                                    view_info_indexer_gcp,
                                    gcps_relative) != 0)
    {
      return -1;
    }

    ExtrinsicParamsContainer extrinsic_params_set_absolute_estimate =
      extrinsic_params_set_relative_estimate;
    size_t number_of_extrinsics =
      extrinsic_params_set_absolute_estimate.size();
    for (size_t i = 0; i < number_of_extrinsics; i++)
    {
      ExtrinsicParams& extrinsic_params =
        extrinsic_params_set_absolute_estimate[i];
      extrinsic_params.rotation() =
        extrinsic_params.rotation() * rotation_similar.Inverse();
      extrinsic_params.position() =
        scale_similar * (rotation_similar * extrinsic_params.position()) +
        translate_similar;
    }

    PointContainer points_absolute_estimate = points_relative_estimate;
    size_t number_of_points = points_absolute_estimate.size();
    for (size_t i = 0; i < number_of_points; i++)
    {
      Point& point = points_absolute_estimate[i];
      point = scale_similar * (rotation_similar * point) + translate_similar;
    }

    PointContainer gcps_absolute_estimate = gcps_relative;
    size_t number_of_available_gcps = gcps_relative.size();
    for (size_t i = 0; i < number_of_available_gcps; i++)
    {
      Point& gcp = gcps_absolute_estimate[i];
      gcp = scale_similar * (rotation_similar * gcp) + translate_similar;
    }

    PointContainer gcps_absolute_true_reordered(number_of_available_gcps);
    size_t number_of_gcps = gcps.size();
    for (size_t i = 0; i < number_of_gcps; i++)
    {
      if (track_point_map_gcp.IsValid(i))
      {
        gcps_absolute_true_reordered[track_point_map_gcp[i]] = gcps[i];
      }
    }

    std::string extrinsic_accuracy_incremental_path =
      test_name_ + "_extrinsic_accuracy_incremental.txt";
    int result = 0;
    std::string gcp_accuracy_incremental_path =
      test_name_ + "_gcp_accuracy_incremental.txt";
    if (tester.TestPointsAccuracy(
        gcps_absolute_true_reordered,
        gcps_absolute_estimate,
        gcp_accuracy_incremental_path,
        ground_resolution_ * 4) != 0) result = -1;

    if (!extrinsic_params_set_absolute_true.empty())
    {

      ExtrinsicParamsContainer extrinsic_params_set_absolute_true_reordered(
        extrinsic_params_set_absolute_estimate.size());
      size_t number_of_images = extrinsic_params_set_absolute_true.size();
      for (size_t i = 0; i < number_of_images; i++)
      {
        if (image_extrinsic_map.IsValid(i))
        {
          extrinsic_params_set_absolute_true_reordered[image_extrinsic_map[i]] =
            extrinsic_params_set_absolute_true[i];
        }
      }

      if (tester.TestExtrinsicAccuracy(
            extrinsic_params_set_absolute_true_reordered,
            extrinsic_params_set_absolute_estimate,
            extrinsic_accuracy_incremental_path,
            ground_resolution_ * 8) != 0) result = -1;
    }

    if (!points_absolute_true.empty())
    {
      PointContainer points_absolute_true_reordered(
                       points_absolute_estimate.size());
      size_t number_of_tracks = points_absolute_true.size();
      for (size_t i = 0; i < number_of_tracks; i++)
      {
        if (track_point_map.IsValid(i))
        {
          points_absolute_true_reordered[track_point_map[i]] =
            points_absolute_true[i];
        }
      }

      std::string point_accuracy_incremental_path =
        test_name_ + "_point_accuracy_incremental.txt";
      if (tester.TestPointsAccuracy(
            points_absolute_true_reordered,
            points_absolute_estimate,
            point_accuracy_incremental_path,
            ground_resolution_ * 4) != 0) result = -1;
    }

    return result;
  }

private:
  Scalar key_stddev_;
  std::string test_name_;
  Scalar ground_resolution_;
};

template <typename _Scalar>
struct RichTrack
{
  typedef _Scalar Scalar;
  hs::sfm::Track track;
  EIGEN_VECTOR(Scalar, 3) point;

  bool operator < (const RichTrack<Scalar>& other) const
  {
    return (track < other.track);
  }
};

TEST(TestIncremental, SyntheticTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestIncremental<Scalar> Tester;
  typedef Tester::IntrinsicParams IntrinsicParams;
  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Tester::ExtrinsicParams ExtrinsicParams;
  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef Tester::KeysetContainer KeysetContainer;
  typedef Tester::PointContainer PointContainer;
  typedef Tester::MatchContainer MatchContainer;
  typedef Tester::TrackContainer TrackContainer;

  typedef RichTrack<Scalar> RichTrack;
  typedef EIGEN_STD_VECTOR(RichTrack) RichTrackContainer;

  typedef hs::sfm::incremental::SyntheticDataGenerator<Scalar, ImageDimension>
          Generator;

  typedef Generator::FlightGenerator FlightGenerator;
  typedef Generator::FlightGeneratorContainer FlightGeneratorContainer;
  typedef Generator::Image Image;
  typedef Generator::ImageContainer ImageContainer;

  FlightGeneratorContainer flight_generators;
  IntrinsicParamsContainer intrinsic_params_set_true;
  IntrinsicParamsContainer intrinsic_params_set_initial;

  Scalar focal_length_in_metre_0 = 0.018858358970276164;
  size_t number_of_strips_0 = 3;
  size_t number_of_cameras_in_strip_0 = 4;
  Scalar ground_resolution_0 = 0.1;
  ImageDimension image_width_0 = 6000;
  ImageDimension image_height_0 = 4000;
  Scalar pixel_size_0 = 0.0000039;
  size_t number_of_points_0 = 10000;
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
                                     3000-42.4095312016,
                                     2000+31.699212823,
                                     1,
                                     -0.0050490462006048831,
                                     0.031293804298609881,
                                     -0.030794820960442223,
                                     -0.00055376548320189127,
                                     -0.00049877717768381476);
  intrinsic_params_set_true.push_back(intrinsic_params_0);
  intrinsic_params_set_initial.push_back(
    IntrinsicParams(4871.79487179487179,
                    0,
                    Scalar(image_width_0 / 2),
                    Scalar(image_height_0 / 2)));

  Scalar focal_length_in_metre_1 = 0.02995452167701055;
  size_t number_of_strips_1 = 3;
  size_t number_of_cameras_in_strip_1 = 4;
  Scalar ground_resolution_1 = 0.1;
  ImageDimension image_width_1 = 6000;
  ImageDimension image_height_1 = 4000;
  Scalar pixel_size_1 = 0.0000039;
  size_t number_of_points_1 = 10000;
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
                                     3000-21.669436058,
                                     2000+44.8644764322,
                                     1,
                                     -0.02529179096221609,
                                     0.23762413973445157,
                                     -0.64208397668697237,
                                     -0.0020605099808780948,
                                     -0.00028706423764766859);
  intrinsic_params_set_true.push_back(intrinsic_params_1);
  intrinsic_params_set_initial.push_back(
    IntrinsicParams(7692.30769230769231,
                    0,
                    Scalar(image_width_1 / 2),
                    Scalar(image_height_1 / 2)));

  Scalar focal_length_in_metre_2 = 0.019056097774998712;
  size_t number_of_strips_2 = 3;
  size_t number_of_cameras_in_strip_2 = 4;
  Scalar ground_resolution_2 = 0.1;
  ImageDimension image_width_2 = 6000;
  ImageDimension image_height_2 = 4000;
  Scalar pixel_size_2 = 0.0000039;
  size_t number_of_points_2 = 10000;
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
                                     3000-35.2052431556,
                                     2000+16.4262220759,
                                     1,
                                     -0.10316088386868619,
                                     0.13520490482776426,
                                     -0.05489235547426094,
                                     4.1434720317373253e-006,
                                     -0.00025018439997095336);
  intrinsic_params_set_true.push_back(intrinsic_params_2);
  intrinsic_params_set_initial.push_back(
    IntrinsicParams(4871.79487179487179,
                    0,
                    Scalar(image_width_2 / 2),
                    Scalar(image_height_2 / 2)));

  Scalar flight_longitudinal_overlap_ratio = 0.9;
  Scalar flight_lateral_overlap_ratio = 0.2;
  Scalar north_west_angle = 0.0;
  Scalar north_west_angle_stddev = 5.0;
  Scalar offset_stddev = 15.0;
  size_t number_of_points = 40000;
  size_t number_of_gcps = 10;
  //TODO:Outlier test needed!
  Scalar outlier_ratio = 0.00;
  Scalar key_stddev = 1.0;

  Generator generator(flight_longitudinal_overlap_ratio,
                      flight_lateral_overlap_ratio,
                      north_west_angle,
                      north_west_angle_stddev,
                      offset_stddev,
                      flight_generators,
                      outlier_ratio,
                      key_stddev,
                      intrinsic_params_set_true,
                      number_of_points);

  ExtrinsicParamsContainer extrinsic_params_set_absolute_true;
  ImageContainer images;
  PointContainer points_absolute;
  KeysetContainer keysets_true;
  TrackContainer tracks;
  hs::sfm::CameraViewContainer camera_views;
  std::vector<size_t> image_intrinsic_map;

  ASSERT_EQ(0, generator.GenerateAbsoluteScene(
                           extrinsic_params_set_absolute_true,
                           images,
                           points_absolute,
                           keysets_true,
                           image_intrinsic_map,
                           tracks,
                           camera_views));
  RichTrackContainer rich_tracks;
  for (size_t i = 0; i < number_of_points; i++)
  {
    if (tracks[i].size() > 1)
    {
      RichTrack rich_track;
      rich_track.point = points_absolute[i];
      rich_track.track = tracks[i];
      std::sort(rich_track.track.begin(), rich_track.track.end());
      rich_tracks.push_back(rich_track);
    }
  }
  std::sort(rich_tracks.begin(), rich_tracks.end());
  size_t number_of_rich_tracks = rich_tracks.size();
  TrackContainer tracks_true(number_of_rich_tracks);
  PointContainer points_absolute_true(number_of_rich_tracks);
  for (size_t i = 0; i < number_of_rich_tracks; i++)
  {
    tracks_true[i] = rich_tracks[i].track;
    points_absolute_true[i] = rich_tracks[i].point;
  }

  hs::sfm::MatchesTracksConvertor matches_tracks_convertor;
  MatchContainer matches;
  ASSERT_EQ(0, matches_tracks_convertor(tracks_true, matches));

  KeysetContainer keysets_noised;
  ASSERT_EQ(0, generator.GenerateNoisedKeysets(keysets_true,
                                               images,
                                               keysets_noised));

  KeysetContainer keysets_gcp_true;
  PointContainer gcps;
  TrackContainer tracks_gcp;
  ASSERT_EQ(0, generator.GenerateGCPData(intrinsic_params_set_true,
                                         extrinsic_params_set_absolute_true,
                                         images,
                                         image_intrinsic_map,
                                         number_of_gcps,
                                         keysets_gcp_true,
                                         gcps,
                                         tracks_gcp));

  KeysetContainer keysets_gcp_noised;
  ASSERT_EQ(0, generator.GenerateNoisedKeysets(keysets_gcp_true,
                                               images,
                                               keysets_gcp_noised));

  std::string test_name = "synthetic_data";

  Tester tester(0.1, key_stddev, test_name);
  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
                      image_intrinsic_map,
                      matches,
                      keysets_noised,
                      gcps,
                      tracks_gcp,
                      keysets_gcp_noised,
                      intrinsic_params_set_true,
                      extrinsic_params_set_absolute_true,
                      points_absolute_true));

}

TEST(TestIncremental, RealTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestIncremental<Scalar> Tester;
  typedef Tester::IntrinsicParams IntrinsicParams;
  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Tester::KeysetContainer KeysetContainer;
  typedef Tester::PointContainer PointContainer;

  typedef hs::sfm::incremental::RealDataGenerator<Scalar> Generator;

  typedef hs::sfm::MatchContainer MatchContainer;
  typedef hs::sfm::TrackContainer TrackContainer;

  std::string data_path = hs::test::getTestDataPath();
  std::string out_path = data_path + "sfm/incremental/real_data/bundler.out";
  std::string gcp_path = data_path + "sfm/incremental/real_data/gcp.xml";
  IntrinsicParamsContainer intrinsic_params_set_initial;
  intrinsic_params_set_initial.push_back(IntrinsicParams(7692.30769230769231,
                                                         0,
                                                         3000,
                                                         2000));
  intrinsic_params_set_initial.push_back(IntrinsicParams(7692.30769230769231,
                                                         0,
                                                         3000,
                                                         2000));
  intrinsic_params_set_initial.push_back(IntrinsicParams(7692.30769230769231,
                                                         0,
                                                         3000,
                                                         2000));
  std::vector<size_t> image_intrinsic_map(123);
  for (size_t i = 0; i < 2; i++)
  {
    image_intrinsic_map[i] = 0;
  }
  for (size_t i = 2; i < 20; i++)
  {
    image_intrinsic_map[i] = 1;
  }
  for (size_t i = 20; i < 123; i++)
  {
    image_intrinsic_map[i] = 2;
  }

  KeysetContainer keysets;
  MatchContainer matches;
  PointContainer gcps;
  TrackContainer tracks_gcp;
  KeysetContainer keysets_gcp;
  ASSERT_EQ(0, Generator::LoadBundlerOutFile(out_path, 6000, 4000,
                                             keysets, matches));
  ASSERT_EQ(0, Generator::LoadGCPs(gcp_path, gcps, tracks_gcp, keysets_gcp));

  Scalar key_stddev = Scalar(1);
  std::string test_name = "real_data";
  Tester tester(0.1, key_stddev, test_name);
  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
                      image_intrinsic_map,
                      matches,
                      keysets,
                      gcps,
                      tracks_gcp,
                      keysets_gcp));
}

}
