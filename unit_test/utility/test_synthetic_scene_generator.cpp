#include <gtest/gtest.h>

#include "hs_sfm/utility/synthetic_scene_generator.hpp"
#include "hs_sfm/utility/sfm_file_io.hpp"

namespace
{

template <typename _Scalar>
struct PtsPrjCamResidual
{
  typedef _Scalar Scalar;
  typedef hs::sfm::IntrinParam<Scalar> Intrin;
  typedef EIGEN_VECTOR(Intrin) IntrinContainer;
  typedef hs::sfm::ExtrinParam<Scalar> Extrin;
  typedef EIGEN_VECTOR(Extrin) ExtrinContainer;
  typedef hs::sfm::ImageKeys<Scalar> ImgKeys;
  typedef EIGEN_VECTOR(ImgKeys) ImgKeysContainer;
  typedef EIGEN_VEC(Scalar, 3) Pt3D;
  typedef EIGEN_VEC(Scalar, 2) Key;
  typedef Pt3D KeyH;
  typedef EIGEN_VECTOR(Pt3D) Pt3DContainer;
  typedef std::vector<std::pair<size_t, size_t> > Track;
  typedef std::vector<Track> TrackContainer;
  typedef hs::sfm::CamFunc<Scalar> Cam;
  typedef typename Cam::PMat PMat;

  Scalar operator()(const IntrinContainer& intrins,
            const ExtrinContainer& extrins,
            const Pt3DContainer& pts,
            const ImgKeysContainer& imgKeysSet,
            const TrackContainer& tracks) const
  {
    size_t ptsNum = pts.size();
    size_t camNum = intrins.size();
    if (intrins.size() != camNum || imgKeysSet.size() != camNum ||
      tracks.size() != ptsNum)
    {
      return Scalar(-1);
    }

    Scalar sumDist = Scalar(0);
    size_t measureNum = 0;
    for (size_t i = 0; i < ptsNum; i++)
    {
      size_t viewSize = tracks[i].size();
      measureNum += viewSize;
      for (size_t j = 0; j < viewSize; j++)
      {
        size_t imgId = tracks[i][j].first;
        size_t keyId = tracks[i][j].second;

        Key key = imgKeysSet[imgId][keyId];

        PMat P = Cam::getPMat(intrins[imgId], extrins[imgId]);
        KeyH keyH = P.block(0, 0, 3, 3) * pts[i] + 
               P.block(0, 3, 3, 1);
        keyH /= keyH(2);

        Scalar dist = (keyH.segment(0, 2) - key).squaredNorm();
        sumDist += dist;
      }
    }

    Scalar resd = sqrt(sumDist / Scalar(measureNum) / Scalar(2));

    return resd;
  }
};

TEST(TestSyntheticSceneGenerator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImgDim;

  typedef hs::sfm::SceneGenerator<Scalar, ImgDim> SceneGen;
  typedef hs::sfm::KeysGenerator<Scalar, ImgDim> KeysGen;
  typedef hs::sfm::MatchGenerator MatchesGen;
  typedef hs::sfm::NoiseImageKeys<Scalar> NoiseImgKeys;
  typedef hs::sfm::NoiseExtrin<Scalar> NoiseExtr;

  typedef SceneGen::Intrin Intrin;
  typedef SceneGen::Extrin Extrin;
  typedef SceneGen::IntrinContainer IntrinContainer;
  typedef SceneGen::ExtrinContainer ExtrinContainer;
  typedef SceneGen::ImageContainer ImageContainer;
  typedef SceneGen::Pt3DContainer Pt3DContainer;

  typedef KeysGen::ImgKeys ImgKeys;
  typedef KeysGen::ImgKeysContainer ImgKeysContainer;
  typedef KeysGen::TrackContainer TrackContainer;

  typedef hs::sfm::MatchContainer MatchContainer;

  typedef NoiseImgKeys::Cov KeyCov;

  typedef NoiseExtr::Cov ExtrinCov;

  typedef std::pair<size_t, size_t> CamPairKey;
  typedef EIGEN_MAP(CamPairKey, Extrin) CamPairContainer;

  typedef PtsPrjCamResidual<Scalar> Residual;

  Scalar f = 0.006;
  size_t stripNum = 15;
  size_t camsNumInStrip = 20;
  Scalar grdRes = 0.1;
  ImgDim imgW = 3648;
  ImgDim imgH = 2736;
  Scalar pixSize = 0.00000203311408298266;
  size_t ptsNum = 2000;
  Scalar lateralOverlap = 0.6;
  Scalar longitudinalOverlap = 0.8;
  Scalar sceneMaxHeight = 50;
  Scalar camHeightDev = 2;
  Scalar camPlannarDev = 2;
  Scalar camRotDev = 10;
  Scalar nwAngle = 60;

  SceneGen sceneGen(f, stripNum, camsNumInStrip, grdRes, 
            imgW, imgH, pixSize, ptsNum,
            lateralOverlap, longitudinalOverlap,
            sceneMaxHeight,
            camHeightDev, camPlannarDev, camRotDev, nwAngle);

  IntrinContainer intrins;
  ExtrinContainer extrins;
  ImageContainer images;
  Pt3DContainer pts;
  sceneGen(intrins, extrins, images, pts);
  //保存文件方便可视化观察
  //sfm::SaveXugFile<Scalar, ImgDim> saver;
  //saver("TestSyntheticSceneGenerator/SimpleTest/test.xug",
  //    intrins, extrins, images, pts, 10);

  KeysGen keysGen(imgW, imgH);
  ImgKeysContainer imgKeysSet;
  TrackContainer tracks;
  keysGen(intrins, extrins, pts, imgKeysSet, tracks);

  MatchesGen mathesGen;
  MatchContainer matches;
  mathesGen(tracks, matches);

  //特征点加入噪声
  ImgKeysContainer noisedImgKeysSet = imgKeysSet;
  size_t numCam = noisedImgKeysSet.size();
  KeyCov keyCov = KeyCov::Identity();
  NoiseImgKeys noiseImgKeys(keyCov);
  for (size_t i = 0; i < numCam; i++)
  {
    noiseImgKeys(noisedImgKeysSet[i]);
  }

  //外参数加入噪声，作为POS信息
  ExtrinContainer priorExtrins = extrins;
  ExtrinCov covR = ExtrinCov::Identity();
  covR *= Scalar(5) / 180 * Scalar(M_PI);
  ExtrinCov covC = ExtrinCov::Identity();
  covC *= Scalar(10);
  NoiseExtr noiseExtr(covR, covC);
  for (size_t i = 0; i < numCam; i++)
  {
    noiseExtr(priorExtrins[i]);
  }

  //保存文件方便可视化观察
  //saver("TestSyntheticSceneGenerator/SimpleTest/test_noise.xug",
  //  intrins, priorExtrins, images, pts, 10);

}

}
