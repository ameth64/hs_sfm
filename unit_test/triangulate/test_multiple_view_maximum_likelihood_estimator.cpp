#include <iostream>
#include <cmath>

#include <gtest/gtest.h>

#include "hs_test_utility/test_env/data_path.hpp"

#include "hs_sfm/triangulate/multiple_view_maximum_likelihood_estimator.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/key_type.hpp"

#include "test_multiple_view_maximum_likelihood_estimator_real_data_generator.hpp"

namespace
{

template <typename _Scalar>
class TestMultipleViewMaximumLikelihoodEstimator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::triangulate::MultipleViewMaximumLikelihoodEstimator<Scalar>
          Estimator;
  typedef typename Estimator::ExtrinsicParams ExtrinsicParams;
  typedef typename Estimator::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef typename Estimator::IntrinsicParams IntrinsicParams;
  typedef typename Estimator::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef typename Estimator::Key Key;
  typedef typename Estimator::KeyContainer KeyContainer;
  typedef typename Estimator::Point Point;

  typedef EIGEN_STD_VECTOR(Point) PointContainer;
  typedef hs::sfm::ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;

public:
  Err operator() (const IntrinsicParamsContainer& intrinsic_params_set,
                  const ExtrinsicParamsContainer& extrinsic_params_set,
                  const KeysetContainer& keysets,
                  const hs::sfm::TrackContainer& tracks,
                  const PointContainer& points_true,
                  const hs::sfm::ObjectIndexMap& image_intrinsic_map,
                  const hs::sfm::ObjectIndexMap& image_extrinsic_map,
                  Scalar mean_planar_threshold,
                  Scalar mean_height_threshold) const
  {
    Scalar mean_planar_error = Scalar(0);
    Scalar mean_height_error = Scalar(0);
    size_t number_of_points = tracks.size();
    if (points_true.size() != number_of_points) return -1;
    size_t number_of_images = image_intrinsic_map.Size();
    if (image_extrinsic_map.Size() != number_of_images) return -1;
    Estimator estimator;
    for (size_t point_id = 0; point_id < number_of_points; point_id++)
    {
      size_t number_of_views = tracks[point_id].size();
      if (number_of_views < 2)
      {
        return -1;
      }
      IntrinsicParamsContainer intrinsic_params_set_triangulate;
      ExtrinsicParamsContainer extrinsic_params_set_triangulate;
      KeyContainer keys_triangulate;
      for (size_t view_id = 0; view_id < number_of_views; view_id++)
      {
        size_t image_id = tracks[point_id][view_id].first;
        size_t key_id = tracks[point_id][view_id].second;
        if (!image_intrinsic_map.IsValid(image_id)) continue;
        if (!image_extrinsic_map.IsValid(image_id)) continue;
        size_t intrinsic_id = image_intrinsic_map[image_id];
        size_t extrinsic_id = image_extrinsic_map[image_id];
        intrinsic_params_set_triangulate.push_back(
          intrinsic_params_set[intrinsic_id]);
        extrinsic_params_set_triangulate.push_back(
          extrinsic_params_set[extrinsic_id]);
        keys_triangulate.push_back(keysets[image_id][key_id]);
      }
      if (intrinsic_params_set_triangulate.size() < 2) return -1;
      Point point_estimate;
      if (estimator(intrinsic_params_set_triangulate,
                    extrinsic_params_set_triangulate,
                    keys_triangulate,
                    point_estimate) != 0) return -1;
      Scalar planar_error = 
        (point_estimate - points_true[point_id]).segment(0, 2).norm();
      Scalar height_error =
        std::abs(point_estimate[2] - points_true[point_id][2]);
      mean_planar_error += planar_error;
      mean_height_error += height_error;
        
    }
    mean_planar_error /= Scalar(number_of_points);
    mean_height_error /= Scalar(number_of_points);

    std::cout<<"mean_planar_error:"<<mean_planar_error<<"\n";
    if (mean_planar_error > mean_planar_threshold)
    {
      return -1;
    }

    std::cout<<"mean_height_error:"<<mean_height_error<<"\n";
    if (mean_height_error > mean_height_threshold)
    {
      return -1;
    }

    return 0;
  }
};

TEST(TestMultipleViewMaximumLikelihoodEstimator, Real242ImagesTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestMultipleViewMaximumLikelihoodEstimator<Scalar> Tester;
  typedef Tester::IntrinsicParams IntrinsicParams;
  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Tester::ExtrinsicParams ExtrinsicParams;
  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef Tester::KeysetContainer KeysetContainer;
  typedef Tester::PointContainer PointContainer;

  typedef hs::sfm::triangulate::
          TestMultipleViewMaximumLikelihoodEstimatorRealDataGenerator<Scalar>
          Generator;
  typedef hs::sfm::TrackContainer TrackContainer;
  std::string data_path = hs::test::getTestDataPath();
  std::string out_path =
    data_path + "sfm/incremental/real_data_242_images/bundler.out";
  std::string gcp_path =
    data_path + "sfm/incremental/real_data_242_images/gcp.xml";
  IntrinsicParamsContainer intrinsic_params_set_initial;
  intrinsic_params_set_initial.push_back(IntrinsicParams(4666.67,
                                                         0,
                                                         3000,
                                                         2000));
  IntrinsicParamsContainer intrinsic_params_set;
  intrinsic_params_set.push_back(IntrinsicParams(4880.22,
                                                 0,
                                                 3026.11,
                                                 1992.36,
                                                 1,
                                                 -0.100423,
                                                 0.128487,
                                                 -0.0482081,
                                                 -3.17902e-05,
                                                 7.75629e-05));
  hs::sfm::ObjectIndexMap image_intrinsic_map(242);
  for (size_t i = 0; i < 242; i++)
  {
    image_intrinsic_map[i] = 0;
  }

  ExtrinsicParamsContainer extrinsic_params_set;
  hs::sfm::ObjectIndexMap image_extrinsic_map;
  KeysetContainer keysets_noised;
  PointContainer points_true;
  TrackContainer tracks;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar key_stddev = Scalar(1);
  Scalar outlier_ratio = Scalar(0.0);

  Generator genertor(outlier_ratio, key_stddev);
  genertor.Generate(out_path,
                    intrinsic_params_set,
                    image_intrinsic_map,
                    extrinsic_params_set,
                    image_extrinsic_map,
                    image_width,
                    image_height,
                    tracks,
                    keysets_noised,
                    points_true);

  Tester tester;
  ASSERT_EQ(0, tester(intrinsic_params_set,
                      extrinsic_params_set,
                      keysets_noised,
                      tracks,
                      points_true,
                      image_intrinsic_map,
                      image_extrinsic_map,
                      0.2, 0.4));
}

TEST(TestMultipleViewMaximumLikelihoodEstimator, Real1067ImagesTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef TestMultipleViewMaximumLikelihoodEstimator<Scalar> Tester;
  typedef Tester::IntrinsicParams IntrinsicParams;
  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Tester::ExtrinsicParams ExtrinsicParams;
  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef Tester::KeysetContainer KeysetContainer;
  typedef Tester::PointContainer PointContainer;

  typedef hs::sfm::triangulate::
          TestMultipleViewMaximumLikelihoodEstimatorRealDataGenerator<Scalar>
          Generator;
  typedef hs::sfm::TrackContainer TrackContainer;
  std::string data_path = hs::test::getTestDataPath();
  std::string out_path =
    data_path + "sfm/incremental/real_data_1067_images/bundler.out";
  std::string gcp_path =
    data_path + "sfm/incremental/real_data_1067_images/gcp.xml";
  IntrinsicParamsContainer intrinsic_params_set;
  intrinsic_params_set.push_back(IntrinsicParams(7680.93,
                                                      0,
                                                      3005.66,
                                                      2002.71,
                                                      1,
                                                      -0.139239,
                                                      0.264343,
                                                      0.157404));
  intrinsic_params_set.push_back(IntrinsicParams(7687.72,
                                                      0,
                                                      3001.87,
                                                      2010.12,
                                                      1,
                                                      -0.138811,
                                                      0.250986,
                                                      0.20556,
                                                      -0.00040312,
                                                      0.000140978));

  hs::sfm::ObjectIndexMap image_intrinsic_map(1067);
  for (size_t i = 0; i < 491; i++)
  {
    image_intrinsic_map[i] = 0;
  }
  for (size_t i = 491; i < 1067; i++)
  {
    image_intrinsic_map[i] = 1;
  }

  ExtrinsicParamsContainer extrinsic_params_set;
  hs::sfm::ObjectIndexMap image_extrinsic_map;
  KeysetContainer keysets_noised;
  PointContainer points_true;
  TrackContainer tracks;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar key_stddev = Scalar(1);
  Scalar outlier_ratio = Scalar(0.0);

  Generator genertor(outlier_ratio, key_stddev);
  genertor.Generate(out_path,
                    intrinsic_params_set,
                    image_intrinsic_map,
                    extrinsic_params_set,
                    image_extrinsic_map,
                    image_width,
                    image_height,
                    tracks,
                    keysets_noised,
                    points_true);

  Tester tester;
  ASSERT_EQ(0, tester(intrinsic_params_set,
                      extrinsic_params_set,
                      keysets_noised,
                      tracks,
                      points_true,
                      image_intrinsic_map,
                      image_extrinsic_map,
                      0.2, 0.4));
}

}