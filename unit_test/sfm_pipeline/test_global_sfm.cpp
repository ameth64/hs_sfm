#include <gtest/gtest.h>

#include "data_tester.hpp"
#include "sfm_pipeline_tester.hpp"
#include "synthetic_data_generator.hpp"
#include "real_data_generator.hpp"
#include "real_synthetic_data_generator.hpp"
#include "synthetic_1dsfm_data_generator.hpp"

#include "hs_test_utility/test_env/data_path.hpp"

#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"
#include "hs_sfm/sfm_pipeline/global_sfm.hpp"
#include "hs_sfm/sfm_pipeline/gcp_similar_transform_estimator.hpp"
#include "hs_sfm/triangulate/multiple_view_maximum_likelihood_estimator.hpp"
#include "hs_sfm/sfm_pipeline/bundle_adjustment_gcp_constrained_optimizor.hpp"
#define DEBUG_TMP 1
#if DEBUG_TMP
#include "hs_sfm/sfm_utility/debug_tmp.hpp"
#include "hs_sfm/sfm_file_io/scene_ply_saver.hpp"
#include <fstream>
#include <iomanip>
#endif

#define TRY_1DSFM 1

namespace
{

TEST(TestGlobalSFM, RealSynthetic242ImageTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;
  typedef hs::sfm::pipeline::GlobalSFM<Scalar> SFMPipeline;
  typedef hs::sfm::pipeline::SFMPipelineTester<SFMPipeline> Tester;
  typedef Tester::IntrinsicParams IntrinsicParams;
  typedef Tester::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Tester::ExtrinsicParams ExtrinsicParams;
  typedef Tester::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef Tester::KeysetContainer KeysetContainer;
  typedef Tester::PointContainer PointContainer;

  typedef hs::sfm::pipeline::RealSyntheticDataGenerator<Scalar> Generator;

  typedef hs::sfm::MatchContainer MatchContainer;
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
  IntrinsicParamsContainer intrinsic_params_set_true;
  //intrinsic_params_set_true.push_back(IntrinsicParams(4880.22,
  //                                                    0,
  //                                                    3026.11,
  //                                                    1992.36,
  //                                                    1,
  //                                                    -0.100423,
  //                                                    0.128487,
  //                                                    -0.0482081,
  //                                                    -3.17902e-05,
  //                                                    7.75629e-05));
  intrinsic_params_set_true.push_back(IntrinsicParams(4666.67,
                                                      0,
                                                      3000,
                                                      2000));
  std::vector<size_t> image_intrinsic_map(242);
  for (size_t i = 0; i < 242; i++)
  {
    image_intrinsic_map[i] = 0;
  }

  KeysetContainer keysets;
  TrackContainer tracks;
  MatchContainer matches;
  PointContainer gcps;
  TrackContainer tracks_gcp;
  KeysetContainer keysets_gcp;
  PointContainer check_points;
  TrackContainer tracks_check_point;
  KeysetContainer keysets_check_point;
  PointContainer points_true;
  ExtrinsicParamsContainer extrinsic_params_set_true;
  Scalar key_stddev = Scalar(1);
  Scalar outlier_ratio = Scalar(0.1);

  Generator generator(outlier_ratio, key_stddev);
  generator.Generate(out_path, gcp_path,
                     intrinsic_params_set_true, image_intrinsic_map,
                     6000, 4000, tracks, matches,
                     points_true, extrinsic_params_set_true, keysets,
                     gcps, tracks_gcp, keysets_gcp,
                     check_points, tracks_check_point, keysets_check_point);

#if TRY_1DSFM
  std::ofstream eg_file("EGs_true.txt");
  eg_file.setf(std::ios::fixed);
  eg_file<<std::setprecision(8);
  auto itr_image_pair = matches.begin();
  auto itr_image_pair_end = matches.end();
  typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  for (; itr_image_pair != itr_image_pair_end; ++itr_image_pair)
  {
    size_t image_id_left = itr_image_pair->first.first;
    size_t image_id_right = itr_image_pair->first.second;
    const ExtrinsicParams& extrinsic_params_left =
      extrinsic_params_set_true[image_id_left];
    const ExtrinsicParams& extrinsic_params_right =
      extrinsic_params_set_true[image_id_right];
    RMatrix R_left = extrinsic_params_left.rotation();
    RMatrix R_right = extrinsic_params_right.rotation();
    Vector3 c_left = extrinsic_params_left.position();
    Vector3 c_right = extrinsic_params_right.position();
    R_left.row(1) *= Scalar(-1);
    R_left.row(2) *= Scalar(-1);
    R_right.row(1) *= Scalar(-1);
    R_right.row(2) *= Scalar(-1);
    RMatrix R = R_left * R_right.transpose();
    Vector3 c = c_right - c_left;
    c /= c.norm();
    c = R_left * c;
    eg_file<<image_id_left<<" "<<image_id_right<<" "
           <<R(0, 0)<<" "<<R(0, 1)<<" "<<R(0, 2)<<" "
           <<R(1, 0)<<" "<<R(1, 1)<<" "<<R(1, 2)<<" "
           <<R(2, 0)<<" "<<R(2, 1)<<" "<<R(2, 2)<<" "
           <<c(0)<<" "<<c(1)<<" "<<c(2)<<"\n";
  }
  eg_file.close();

  //Output gt_bundle.out
  std::ofstream gt_bundle_file("gt_bundle.out");
  gt_bundle_file.setf(std::ios::fixed);
  gt_bundle_file<<std::setprecision(10);
  gt_bundle_file<<"# Bundle file v0.3\n";
  gt_bundle_file<<extrinsic_params_set_true.size()<<" "
                <<points_true.size()<<"\n";
  RMatrix R0 = extrinsic_params_set_true[0].rotation();
  R0.row(1) *= Scalar(-1);
  R0.row(2) *= Scalar(-1);
  Vector3 c0 = extrinsic_params_set_true[0].position();
  for (size_t i = 0; i < extrinsic_params_set_true.size(); i++)
  {
    size_t intrinsic_id = image_intrinsic_map[i];
    const IntrinsicParams& intrinsic_params =
      intrinsic_params_set_true[intrinsic_id];
    const ExtrinsicParams& extrinsic_params = extrinsic_params_set_true[i];
    RMatrix R = extrinsic_params.rotation();
    Vector3 c = extrinsic_params.position();
    R.row(1) *= Scalar(-1);
    R.row(2) *= Scalar(-1);
    R = R0 * R.transpose();
    R.transposeInPlace();
    c = c - c0;
    c = R0 * c;
    Vector3 t = -R * c;
    gt_bundle_file<<intrinsic_params.focal_length()<<" "
                  <<intrinsic_params.k1()<<" "
                  <<intrinsic_params.k2()<<"\n";
    gt_bundle_file<<R(0, 0)<<" "<<R(0, 1)<<" "<<R(0, 2)<<"\n"
                  <<R(1, 0)<<" "<<R(1, 1)<<" "<<R(1, 2)<<"\n"
                  <<R(2, 0)<<" "<<R(2, 1)<<" "<<R(2, 2)<<"\n";
    gt_bundle_file<<t(0)<<" "<<t(1)<<" "<<t(2)<<"\n";
  }
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;
  for (size_t i = 0; i < points_true.size(); i++)
  {
    gt_bundle_file<<points_true[i][0]<<" "
                  <<points_true[i][1]<<" "
                  <<points_true[i][2]<<"\n";
    gt_bundle_file<<"255 255 255\n";
    gt_bundle_file<<tracks[i].size();
    for (size_t j = 0; j < tracks[i].size(); j++)
    {
      size_t image_id = tracks[i][j].first;
      size_t key_id = tracks[i][j].second;
      gt_bundle_file<<" "<<image_id<<" "<<key_id;
      Vector2 key = keysets[image_id][key_id];
      key[0] = key[0] - Scalar(3000);
      key[1] = Scalar(2000) - key[1];
      gt_bundle_file<<" "<<key[0]<<" "<<key[1];
    }
    gt_bundle_file<<"\n";
  }
  gt_bundle_file.close();
#endif

  std::string test_name = "global_real_synthetic_data_242_images";
  SFMPipeline sfm_pipeline;
  Tester tester(sfm_pipeline, 0.1, key_stddev, test_name);
  ASSERT_EQ(0, tester(intrinsic_params_set_initial,
                      image_intrinsic_map,
                      matches,
                      keysets,
                      gcps,
                      tracks_gcp,
                      keysets_gcp,
                      check_points,
                      tracks_check_point,
                      keysets_check_point,
                      intrinsic_params_set_true,
                      extrinsic_params_set_true,
                      points_true));
}

TEST(TestGlobalSFM, 1DSFMTest)
{
  typedef double Scalar;
  typedef hs::sfm::pipeline::Synthetic1DSFMDataGenerator<Scalar> Generator;
  typedef Generator::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef Generator::KeysetContainer KeysetContainer;
  typedef hs::sfm::pipeline::GlobalSFM<Scalar> GlobalSFM;
  typedef GlobalSFM::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef GlobalSFM::PointContainer PointContainer;

  std::string dataset_dir = "/Users/limbo/workspace/test/try_1dsfm/datasets/Ellis_Island";

  Generator generator;
  IntrinsicParamsContainer intrinsic_params_set;
  hs::sfm::ObjectIndexMap image_intrinsic_map;
  hs::sfm::TrackContainer tracks;
  KeysetContainer keysets;
  generator.Generate(dataset_dir,
                     intrinsic_params_set,
                     image_intrinsic_map,
                     tracks,
                     keysets);

  hs::sfm::MatchesTracksConvertor matches_tracks_convertor;
  hs::sfm::MatchContainer matches;
  matches_tracks_convertor(tracks, matches);

  GlobalSFM global_sfm;
  ExtrinsicParamsContainer extrinsic_params_set;
  hs::sfm::ObjectIndexMap image_extrinsic_map;
  PointContainer points;
  hs::sfm::TrackContainer tracks_bundle;
  hs::sfm::ObjectIndexMap track_point_map;
  hs::sfm::ViewInfoIndexer view_info_indexer;
  global_sfm(image_intrinsic_map, matches, keysets,
             intrinsic_params_set, extrinsic_params_set,
             image_extrinsic_map, points, tracks_bundle, track_point_map,
             view_info_indexer);
}

}
