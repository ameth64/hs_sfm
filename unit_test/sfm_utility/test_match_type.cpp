#include <fstream>

#include <gtest/gtest.h>

#include <cereal/archives/portable_binary.hpp>

#include "hs_test_utility/test_toolkit/test_serialization.hpp"

#include "hs_sfm/sfm_utility/match_type.hpp"

namespace
{

TEST(TestMatchType, SerializationTest)
{
  typedef hs::sfm::ImagePair ImagePair;
  typedef hs::sfm::KeyPair KeyPair;
  typedef hs::sfm::MatchContainer MatchContainer;
  typedef hs::sfm::Track Track;
  typedef hs::sfm::TrackContainer TrackContainer;
  typedef hs::sfm::ObjectIndexMap ObjectIndexMap;
  typedef hs::sfm::ViewInfoIndexer ViewInfoIndexer;

  MatchContainer matches;
  matches[ImagePair(0, 1)].push_back(KeyPair(0, 1));
  matches[ImagePair(0, 1)].push_back(KeyPair(0, 2));
  matches[ImagePair(0, 1)].push_back(KeyPair(0, 3));
  matches[ImagePair(0, 1)].push_back(KeyPair(1, 2));
  matches[ImagePair(0, 1)].push_back(KeyPair(1, 3));
  matches[ImagePair(1, 2)].push_back(KeyPair(0, 1));
  matches[ImagePair(1, 2)].push_back(KeyPair(0, 2));
  matches[ImagePair(1, 2)].push_back(KeyPair(0, 3));
  hs::test::TestLocalSerialization<MatchContainer> matches_tester;
  ASSERT_EQ(true, matches_tester(matches, "matches.bin"));

  TrackContainer tracks(100);
  for (size_t i = 0; i < 100; i++)
  {
    for (size_t j = 0; j < 9; j++)
    {
      tracks[i].push_back(std::pair<size_t, size_t>(2, 3));
    }
  }
  hs::test::TestLocalSerialization<TrackContainer> tracks_tester;
  ASSERT_EQ(true, tracks_tester(tracks, "tracks.bin"));

  ObjectIndexMap object_index_map(100);

  for (size_t i = 0; i < 88; i++)
  {
    object_index_map[i + 12] = i;
  }
  hs::test::TestLocalSerialization<ObjectIndexMap> object_index_map_tester;
  ASSERT_EQ(true, object_index_map_tester(object_index_map,
                                          "object_index_map.bin"));

  ViewInfoIndexer view_info_indexer;
  view_info_indexer.SetViewInfoByTracks(tracks);
  hs::test::TestLocalSerialization<ViewInfoIndexer> view_info_indexer_tester;
  ASSERT_EQ(true, view_info_indexer_tester(view_info_indexer,
                                           "view_info_indexer.bin"));
}

}