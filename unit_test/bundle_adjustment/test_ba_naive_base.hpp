#ifndef _UNIT_TEST_TEST_BA_NAIVE_BASE_HPP_
#define _UNIT_TEST_TEST_BA_NAIVE_BASE_HPP_

#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"
#include "hs_sfm/utility/synthetic_scene_generator.hpp"

template <typename _Scalar, typename _ImgDim>
class SyntheticDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImgDim ImgDim;

  typedef hs::sfm::SceneGenerator<Scalar, ImgDim> SceneGen;
  typedef typename SceneGen::Intrin Intrin;
  typedef typename SceneGen::IntrinContainer IntrinContainer;
  typedef typename SceneGen::Extrin Extrin;
  typedef typename SceneGen::ExtrinContainer ExtrinContainer;
  typedef typename SceneGen::Image Image;
  typedef typename SceneGen::ImageContainer ImageContainer;
  typedef typename SceneGen::Pt3D Pt3D;
  typedef typename SceneGen::Pt3DContainer Pt3DContainer;

  typedef hs::sfm::KeysGenerator<Scalar, ImgDim> KeysGen;
  typedef typename KeysGen::ImgKeys ImgKeys;
  typedef typename KeysGen::ImgKeysContainer ImgKeysContainer;
  typedef typename KeysGen::Track Track;
  typedef typename KeysGen::TrackContainer TrackContainer;

  typedef hs::sfm::ba::BANaiveVecFunc<Scalar> BAVecFunc;
  typedef typename BAVecFunc::Index Index;
  typedef typename BAVecFunc::XVec XVec;
  typedef typename BAVecFunc::YVec YVec;
  typedef typename BAVecFunc::FeatMap FeatMap;
  typedef typename BAVecFunc::FeatMapContainer FeatMapContainer;
  typedef typename BAVecFunc::Vec3 Vec3;

  SyntheticDataGenerator(Scalar f, size_t stripNum,
              size_t camsNumInStrip, Scalar grdRes,
              ImgDim imgW, ImgDim imgH, Scalar pixSize, size_t ptsNum,
              Scalar lateralOverlap,
              Scalar longitudinalOverlap,
              Scalar sceneMaxHeight,
              Scalar camHeightDev,
              Scalar camPlannarDev,
              Scalar camRotDev,
              Scalar nwAngle)
    : m_sceneGen(f, stripNum, camsNumInStrip, grdRes,
            imgW, imgH, pixSize, ptsNum, 
            lateralOverlap, longitudinalOverlap, sceneMaxHeight,
            camHeightDev, camPlannarDev, camRotDev, nwAngle),
      m_keysGen(imgW, imgH) {}

  int operator () (BAVecFunc& baVecFunc, XVec& x, YVec& y) const
  {
    IntrinContainer intrins;
    ExtrinContainer extrins;
    ImageContainer images;
    Pt3DContainer pts;
    if (m_sceneGen(intrins, extrins, images, pts) != 0) 
      return -1;
    Scalar f = m_sceneGen.getFInPix();

    ImgKeysContainer imgKeys;
    TrackContainer tracks;
    if (m_keysGen(intrins, extrins, pts, imgKeys, tracks) != 0) 
      return -1;

    size_t keyNum = 0;
    auto itrImgKeys = imgKeys.begin();
    auto itrImgKeysEnd = imgKeys.end();
    std::vector<Index> imgKeysIdOffsets(extrins.size());
    Index i = 0;
    for (; itrImgKeys != itrImgKeysEnd; ++itrImgKeys, ++i)
    {
      imgKeysIdOffsets[i] = keyNum;
      keyNum += itrImgKeys->size();
    }
    auto itrTrack = tracks.begin();
    auto itrTrackEnd = tracks.end();
    FeatMapContainer featMaps(keyNum);
    i = 0;
    for (; itrTrack != itrTrackEnd; ++itrTrack, ++i)
    {
      auto itrView = itrTrack->begin();
      auto itrViewEnd = itrTrack->end();
      for (; itrView != itrViewEnd; ++itrView)
      {
        Index featId = imgKeysIdOffsets[itrView->first] + 
                  itrView->second;
        featMaps[featId].first = Index(itrView->first);
        featMaps[featId].second = Index(i);
      }
    }

    Index camNum = Index(extrins.size());
    Index ptNum = Index(pts.size());
    Index featNum = Index(keyNum);
    baVecFunc.setCamNum(camNum);
    baVecFunc.setPtNum(ptNum);
    baVecFunc.setFeatNum(featNum);
    baVecFunc.setFeatMaps(featMaps);

    //相机外参数
    Index xSize = baVecFunc.getXSize();
    x.resize(xSize);
    Index ySize = baVecFunc.getYSize();
    y.resize(ySize);

    auto extItr = extrins.begin();
    auto extItrEnd = extrins.end();
    i = 0;
    for (; extItr != extItrEnd; ++extItr, ++i)
    {
      Vec3 t = -(extItr->m_r * extItr->m_c);
      x[i * BAVecFunc::m_paramsPerCam + 0] = extItr->m_r[0];
      x[i * BAVecFunc::m_paramsPerCam + 1] = extItr->m_r[1];
      x[i * BAVecFunc::m_paramsPerCam + 2] = extItr->m_r[2];
      x[i * BAVecFunc::m_paramsPerCam + 3] = t[0];
      x[i * BAVecFunc::m_paramsPerCam + 4] = t[1];
      x[i * BAVecFunc::m_paramsPerCam + 5] = t[2];
    }

    //点云
    Index camParamSize = baVecFunc.getCamParamsSize();
    auto ptItr = pts.begin();
    auto ptItrEnd = pts.end();
    i = 0;
    for (; ptItr != ptItrEnd; ++ptItr, ++i)
    {
      x.segment(camParamSize + i * BAVecFunc::m_paramsPerPt,
            BAVecFunc::m_paramsPerPt) = pts[i];
    }

    //特征点
    itrImgKeys = imgKeys.begin();
    i = 0;
    for (; itrImgKeys != itrImgKeysEnd; ++itrImgKeys)
    {
      for (size_t keyId = 0; keyId < itrImgKeys->size(); keyId++)
      {
        y.segment(i * BAVecFunc::m_paramsPerFeat,
          BAVecFunc::m_paramsPerFeat) =
          (*itrImgKeys)[keyId] / f;
        i++;
      }
    }

    return 0;
  }

private:
  SceneGen m_sceneGen;
  KeysGen m_keysGen;
};

#endif
