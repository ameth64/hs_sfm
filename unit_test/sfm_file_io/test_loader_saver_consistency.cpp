#include <iostream>

#include <gtest/gtest.h>

#include "hs_sfm/sfm_file_io/keyset_loader.hpp"
#include "hs_sfm/sfm_file_io/keyset_saver.hpp"
#include "hs_sfm/sfm_file_io/matches_loader.hpp"
#include "hs_sfm/sfm_file_io/matches_saver.hpp"
#include "hs_sfm/sfm_file_io/tracks_loader.hpp"
#include "hs_sfm/sfm_file_io/tracks_saver.hpp"
#include "hs_sfm/sfm_file_io/intrinsic_params_set_loader.hpp"
#include "hs_sfm/sfm_file_io/intrinsic_params_set_saver.hpp"
#include "hs_sfm/sfm_file_io/extrinsic_params_set_loader.hpp"
#include "hs_sfm/sfm_file_io/extrinsic_params_set_saver.hpp"
#include "hs_sfm/sfm_file_io/object_index_loader.hpp"
#include "hs_sfm/sfm_file_io/object_index_saver.hpp"

namespace
{

template <typename Loader, typename Saver, typename Data>
int TestConsisitency(const std::string& path)
{
  Data data;
  Loader loader;
  Saver saver;
  if (loader(path, data) != 0) return -1;
  std::string path_duplicate(path);
  size_t pos = path_duplicate.rfind('.');
  if (pos == std::string::npos)
  {
    path_duplicate += "_dup";
  }
  else
  {
    path_duplicate.insert(pos, "_dup");
  }
  if (saver(path_duplicate, data) != 0) return -1;
  Data data_duplicate;
  if (loader(path_duplicate, data_duplicate) != 0) return -1;

  if (data == data_duplicate)
  {
    return 0;
  }
  else
  {
    return -1;
  }
}

struct TracksLoaderData
{
  hs::sfm::TrackContainer tracks;
  hs::sfm::ViewInfoIndexer view_info_indexer;

  bool operator == (const TracksLoaderData& other) const
  {
    return (tracks == other.tracks &&
            view_info_indexer == other.view_info_indexer);
  }
};

struct TracksLoaderWrapper
{
  typedef int Err;
  Err operator() (const std::string& tracks_path,
                  TracksLoaderData& tracks_data) const
  {
    hs::sfm::fileio::TracksLoader tracks_loader;
    return tracks_loader(tracks_path,
                         tracks_data.tracks,
                         tracks_data.view_info_indexer);
  }
};

struct TracksSaverWrapper
{
  typedef int Err;
  Err operator() (const std::string& tracks_path,
                  const TracksLoaderData& tracks_data) const
  {
    hs::sfm::fileio::TracksSaver tracks_saver;
    return tracks_saver(tracks_path,
                        tracks_data.tracks,
                        tracks_data.view_info_indexer);
  }
};

template <typename _Scalar>
class TestLoaderSaverConsistency
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::fileio::KeysetLoader<Scalar> KeysetLoader;
  typedef hs::sfm::fileio::KeysetSaver<Scalar> KeysetSaver;
  typedef typename KeysetLoader::Keys Keyset;

  typedef hs::sfm::fileio::MatchesLoader MatchesLoader;
  typedef hs::sfm::fileio::MatchesSaver MatchesSaver;
  typedef MatchesLoader::MatchContainer MatchContainer;

  typedef hs::sfm::fileio::IntrinsicParamsSetLoader<Scalar> IntrinsicSetLoader;
  typedef hs::sfm::fileio::IntrinsicParamsSetSaver<Scalar> IntrinsicSetSaver;
  typedef typename IntrinsicSetLoader::IntrinsicParamsContainer IntrinsicSet;

  typedef hs::sfm::fileio::ExtrinsicParamsSetLoader<Scalar> ExtrinsicSetLoader;
  typedef hs::sfm::fileio::ExtrinsicParamsSetSaver<Scalar> ExtrinsicSetSaver;
  typedef typename ExtrinsicSetLoader::ExtrinsicParamsContainer ExtrinsicSet;

  typedef hs::sfm::fileio::ObjectIndexLoader ObjectIndexLoader;
  typedef hs::sfm::fileio::ObjectIndexSaver ObjectIndexSaver;
  typedef ObjectIndexLoader::ObjectIndexMap ObjectIndexMap;

public:
  TestLoaderSaverConsistency(
    const std::string& keyset_path,
    const std::string& matches_path,
    const std::string& tracks_path,
    const std::string& intrinsic_set_path,
    const std::string& extrinsic_set_path,
    const std::string& image_extrinsic_path,
    const std::string& track_point_path)
    : keyset_path_(keyset_path),
      matches_path_(matches_path),
      tracks_path_(tracks_path),
      intrinsic_set_path_(intrinsic_set_path),
      extrinsic_set_path_(extrinsic_set_path),
      image_extrinsic_path_(image_extrinsic_path),
      track_point_path_(track_point_path) {}

public:
  Err Test() const
  {
    if (TestKeyset() != 0)
    {
      std::cout<<"test keyset failed!\n";
      return -1;
    }
    if (TestMatches() != 0)
    {
      std::cout<<"test matches failed!\n";
      return -1;
    }
    if (TestTracks() != 0)
    {
      std::cout<<"test tracks failed!\n";
      return -1;
    }
    if (TestIntrinsicSet() != 0)
    {
      std::cout<<"test intrinsic set failed!\n";
      return -1;
    }
    if (TestExtrinsicSet() != 0)
    {
      std::cout<<"test extrinsic set failed!\n";
      return -1;
    }
    if (TestObjectIndex() != 0)
    {
      std::cout<<"test object index failed!\n";
      return -1;
    }

    return 0;
  }

private:
  Err TestKeyset() const
  {
    return TestConsisitency<KeysetLoader, KeysetSaver, Keyset>(keyset_path_);
  }

  Err TestMatches() const
  {
    return TestConsisitency<MatchesLoader, MatchesSaver, MatchContainer>(
             matches_path_);
  }

  Err TestTracks() const
  {
    return TestConsisitency<TracksLoaderWrapper, TracksSaverWrapper,
                            TracksLoaderData>(tracks_path_);
  }

  Err TestIntrinsicSet() const
  {
    return TestConsisitency<IntrinsicSetLoader, IntrinsicSetSaver,
                            IntrinsicSet>(intrinsic_set_path_);
  }

  Err TestExtrinsicSet() const
  {
    return TestConsisitency<ExtrinsicSetLoader, ExtrinsicSetSaver,
                            ExtrinsicSet>(extrinsic_set_path_);
  }

  Err TestObjectIndex() const
  {
    if (TestConsisitency<ObjectIndexLoader, ObjectIndexSaver,
                         ObjectIndexMap>(image_extrinsic_path_) != 0)
      return -1;

    if (TestConsisitency<ObjectIndexLoader, ObjectIndexSaver,
                         ObjectIndexMap>(track_point_path_) != 0)
      return -1;

    return 0;
  }

private:
  std::string keyset_path_;
  std::string matches_path_;
  std::string tracks_path_;
  std::string intrinsic_set_path_;
  std::string extrinsic_set_path_;
  std::string image_extrinsic_path_;
  std::string track_point_path_;
};

TEST(TestLoaderSaverConsistency, SimpleTest)
{
  typedef double Scalar;
  typedef TestLoaderSaverConsistency<Scalar> Tester;
  std::string keyset_path =
    "../../test_data/sfm_file_io/TestLoaderSaverConsistency/IMG_127816_DxO.feat";
  std::string matches_path =
    "../../test_data/sfm_file_io/TestLoaderSaverConsistency/matches.init.txt";
  std::string tracks_path =
    "../../test_data/sfm_file_io/TestLoaderSaverConsistency/tracks.txt";
  std::string intrinsic_set_path =
    "../../test_data/sfm_file_io/TestLoaderSaverConsistency/intrinsic_params_set.txt";
  std::string extrinsic_set_path =
    "../../test_data/sfm_file_io/TestLoaderSaverConsistency/extrinsic_params_set.txt";
  std::string image_extrinsic_path =
    "../../test_data/sfm_file_io/TestLoaderSaverConsistency/image_extrinsic_map.txt";
  std::string track_point_path =
    "../../test_data/sfm_file_io/TestLoaderSaverConsistency/track_point_map.txt";
  Tester tester(keyset_path,
                matches_path,
                tracks_path,
                intrinsic_set_path,
                extrinsic_set_path,
                image_extrinsic_path,
                track_point_path);
  ASSERT_EQ(0, tester.Test());
}

}
