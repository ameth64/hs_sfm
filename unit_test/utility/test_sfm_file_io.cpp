#include <gtest/gtest.h>

#include "hs_test_utility/test_env/data_path.hpp"

#include "hs_sfm/utility/sfm_file_io.hpp"

namespace
{

template <typename Loader, typename Saver, typename Data>
int testIO(const std::string& path)
{
  Data data;
  Loader loader;
  Saver saver;
  if (loader(path, data) != 0) return -1;
  std::string pathDup(path);
  size_t pos = pathDup.rfind('.');
  if (pos == std::string::npos)
  {
    pathDup += "_dup";
  }
  else
  {
    pathDup.insert(pos, "_dup");
  }
  if (saver(pathDup, data) != 0) return -1;
  Data dataDup;
  if (loader(pathDup, dataDup) != 0) return -1;

  if (data == dataDup)
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
  typedef LoadMatches::ImgPair ImgPair;
  typedef LoadMatches::KeyPairContainer KeyPairContainer;
  typedef LoadMatches::MatchContainer MatchContainer;

  typedef typename LoadImgKeys::ImgKeys ImgKeys;

  typedef int Err;

  TestSFMFileIO(const std::string& keyPath,
          const std::string& matchesPath) 
    : m_keyPath(keyPath), 
      m_matchesPath(matchesPath) {}

  Err testImgKeysIO() const
  {
    return testIO<LoadImgKeys, SaveImgKeys, ImgKeys>
      (m_keyPath);
  }

  Err testMatchesIO()
  {
    return testIO<LoadMatches, SaveMatches, MatchContainer>
      (m_matchesPath);
  }

  Err test()
  {
    if (testImgKeysIO() != 0) return -1;
    if (testMatchesIO() != 0) return -1;

    return 0;
  }

  std::string m_keyPath;
  std::string m_matchesPath;
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
