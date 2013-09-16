#include <gtest/gtest.h>

#include "hs_test_utility/test_env/data_path.hpp"

#include "hs_sfm/utility/sfm_file_io.hpp"

namespace
{

template <typename Loader, typename Saver, typename Data>
int TestIO(const std::string& path)
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

template <typename _Scalar>
class TestSFMFileIO
{
public:
  typedef _Scalar Scalar;
  typedef hs::sfm::LoadImageKeys<Scalar> LoadImgKeys;
  typedef hs::sfm::SaveImageKeys<Scalar> SaveImgKeys;

  typedef hs::sfm::LoadMatches LoadMatches;
  typedef hs::sfm::SaveMatches SaveMatches;
  typedef LoadMatches::KeyPair KeyPair;
  typedef LoadMatches::ImagePair ImagePair;
  typedef LoadMatches::KeyPairContainer KeyPairContainer;
  typedef LoadMatches::MatchContainer MatchContainer;

  typedef typename LoadImgKeys::Keys ImgKeys;

  typedef int Err;

  TestSFMFileIO(const std::string& key_path,
                const std::string& matches_path) 
    : key_path_(key_path), 
      matches_path_(matches_path) {}

  Err TestImgKeysIO() const
  {
    return TestIO<LoadImgKeys, SaveImgKeys, ImgKeys>
      (key_path_);
  }

  Err TestMatchesIO()
  {
    return TestIO<LoadMatches, SaveMatches, MatchContainer>
      (matches_path_);
  }

  Err test()
  {
    if (TestImgKeysIO() != 0) return -1;
    if (TestMatchesIO() != 0) return -1;

    return 0;
  }

  std::string key_path_;
  std::string matches_path_;
};

TEST(TestSFMFileIO, SimpleTest)
{
  typedef double Scalar;
  typedef TestSFMFileIO<Scalar> Test;
  std::string dataPath = hs::test::getTestDataPath();
  std::string keyPath = dataPath + 
    "sfm/utility/TestSFMFileIO/IMG_127816_DxO.feat";
  std::string matchesPath = dataPath + 
    "sfm/utility/TestSFMFileIO/matches.init.txt";
  Test test(keyPath, matchesPath);
  ASSERT_EQ(0, test.test());
}

}
