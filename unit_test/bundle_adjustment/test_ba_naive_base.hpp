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

  typedef int Err;

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
  typedef typename KeysGen::CamViewContainer CamViewContainer;

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

  Err operator () (BAVecFunc& baVecFunc, XVec& x, YVec& y) const
  {
    IntrinContainer intrins;
    ExtrinContainer extrins;
    ImageContainer images;
    Pt3DContainer pts;
    if (m_sceneGen(intrins, extrins, images, pts) != 0) 
      return -1;
    Scalar f = getFInPix();

    ImgKeysContainer imgKeys;
    TrackContainer tracks;
    CamViewContainer camViews;
    if (m_keysGen(intrins, extrins, pts, imgKeys, tracks, camViews) != 0) 
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

    std::vector<size_t> camMap(extrins.size(), 0);
    auto itrCamView = camViews.begin();
    auto itrCamViewEnd = camViews.end();
    i = 0;
    size_t viewKeyCamNum = 0;
    for (; itrCamView != itrCamViewEnd; ++itrCamView, ++i)
    {
      if (!itrCamView->empty())
      {
        camMap[i] = viewKeyCamNum + 1;
        viewKeyCamNum++;
      }
    }

    std::vector<size_t> ptMap(pts.size(), 0);
    auto itrTrack = tracks.begin();
    auto itrTrackEnd = tracks.end();
    i = 0;
    size_t viewKeyPtNum = 0;
    for (; itrTrack != itrTrackEnd; ++itrTrack, ++i)
    {
      if (!itrTrack->empty())
      {
        ptMap[i] = viewKeyPtNum + 1;
        viewKeyPtNum++;
      }
    }

    itrTrack = tracks.begin();
    itrTrackEnd = tracks.end();
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
        featMaps[featId].first = Index(camMap[itrView->first] - 1);
        featMaps[featId].second = Index(ptMap[i] - 1);
      }
    }

    Index camNum = Index(viewKeyCamNum);
    Index ptNum = Index(viewKeyPtNum);
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
    itrCamView = camViews.begin();
    itrCamViewEnd = camViews.end();
    i = 0;
    for (; extItr != extItrEnd; ++extItr, ++itrCamView)
    {
      if (!itrCamView->empty())
      {
        Vec3 t = -(extItr->m_r * extItr->m_c);
        x[i * BAVecFunc::m_paramsPerCam + 0] = extItr->m_r[0];
        x[i * BAVecFunc::m_paramsPerCam + 1] = extItr->m_r[1];
        x[i * BAVecFunc::m_paramsPerCam + 2] = extItr->m_r[2];
        x[i * BAVecFunc::m_paramsPerCam + 3] = t[0];
        x[i * BAVecFunc::m_paramsPerCam + 4] = t[1];
        x[i * BAVecFunc::m_paramsPerCam + 5] = t[2];
        i++;
      }
    }

    //点云
    Index camParamSize = baVecFunc.getCamParamsSize();
    auto ptItr = pts.begin();
    auto ptItrEnd = pts.end();
    itrTrack = tracks.begin();
    itrTrackEnd = tracks.end();
    i = 0;
    for (; ptItr != ptItrEnd; ++ptItr, itrTrack++)
    {
      if (!itrTrack->empty())
      {
        x.segment(camParamSize + i * BAVecFunc::m_paramsPerPt,
          BAVecFunc::m_paramsPerPt) = *ptItr;
        i++;
      }
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

  inline Scalar getFInPix() const
  {
    return m_sceneGen.getFInPix();
  }

  Err genGCPs(size_t gcpNum, const BAVecFunc& baVecFunc, const XVec& x,
              Pt3DContainer& gcps,
              ImgKeysContainer& imgKeysSet,
              TrackContainer& tracks) const
  {
    m_sceneGen.genScenePts(gcpNum, gcps);
    Index camNum = baVecFunc.getCamNum();

    IntrinContainer intrins(camNum, Intrin(getFInPix()));
    ExtrinContainer extrins(camNum);

    for (Index i = 0; i < camNum; i++)
    {
      extrins[i].m_r[0] = x[i * BAVecFunc::m_paramsPerCam + 0];
      extrins[i].m_r[1] = x[i * BAVecFunc::m_paramsPerCam + 1];
      extrins[i].m_r[2] = x[i * BAVecFunc::m_paramsPerCam + 2];
      Vec3 t = x.segment(i * BAVecFunc::m_paramsPerCam + 3, 3);
      extrins[i].m_c = -(extrins[i].m_r.inverse() * t);
    }

    CamViewContainer camViews;
    return m_keysGen(intrins, extrins, gcps, imgKeysSet, tracks, camViews);
  }

private:
  SceneGen m_sceneGen;
  KeysGen m_keysGen;
};

#endif
