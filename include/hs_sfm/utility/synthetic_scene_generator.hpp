#ifndef _HS_SFM_UTILITY_SYNTHETIC_SCENE_GENERATOR_HPP_
#define _HS_SFM_UTILITY_SYNTHETIC_SCENE_GENERATOR_HPP_

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_math/geometry/euler_angles.hpp"

#include "hs_sfm/utility/cam_type.hpp"
#include "hs_sfm/utility/key_type.hpp"
#include "hs_sfm/utility/match_type.hpp"
#include "hs_sfm/utility/sfm_file_io.hpp"

namespace hs
{
namespace sfm
{

template <typename _Scalar, typename _ImgDim>
class SceneGenerator
{
public:
typedef _Scalar Scalar;
typedef _ImgDim ImgDim;
typedef IntrinParam<Scalar> Intrin;
typedef EIGEN_VECTOR(Intrin) IntrinContainer;
typedef ExtrinParam<Scalar> Extrin;
typedef typename Extrin::Pos Pos;
typedef typename Extrin::Rot Rot;
typedef EIGEN_VECTOR(Extrin) ExtrinContainer;
typedef ImageParam<ImgDim> Image;
typedef EIGEN_VECTOR(Image) ImageContainer;
typedef EIGEN_VEC(Scalar, 3) Pt3D;
typedef EIGEN_VECTOR(Pt3D) Pt3DContainer;
typedef EIGEN_MAT(Scalar, 3, 3) Mat33;

typedef int Err;

SceneGenerator(Scalar f, size_t stripNum,
     size_t camsNumInStrip, Scalar grdRes,
     ImgDim imgW, ImgDim imgH, Scalar pixSize, size_t ptsNum,
     Scalar lateralOverlap,
     Scalar longitudinalOverlap,
     Scalar sceneMaxHeight,
     Scalar camHeightDev,
     Scalar camPlannarDev,
     Scalar camRotDev,
     Scalar nwAngle)
     : m_f(f), m_stripsNum(stripNum),
     m_camsNumInStrip(camsNumInStrip), 
     m_grdRes(grdRes), m_imgW(imgW), m_imgH(imgH),
     m_pixSize(pixSize),
     m_ptsNum(ptsNum),
     m_lateralOverlap(lateralOverlap),
     m_longitudinalOverlap(longitudinalOverlap),
     m_sceneMaxHeight(sceneMaxHeight),
     m_camHeightDev(camHeightDev),
     m_camPlannarDev(camPlannarDev),
     m_camRotDev(camRotDev),
     m_nwAngle(nwAngle) {}

Scalar getFInPix() const {return m_f / m_pixSize;}

Err operator ()(IntrinContainer& intrins, ExtrinContainer& extrins,
    ImageContainer& images,
    Pt3DContainer& pts) const
{
  Scalar flightHeight = m_f * m_grdRes / m_pixSize;
  Scalar longitudinalRangePerCam =
  (1 - m_longitudinalOverlap) * m_grdRes * Scalar(m_imgH);
  Scalar lateralRangePerCam = 
  (1 - m_lateralOverlap) * m_grdRes * Scalar(m_imgW);

  size_t camNum = m_stripsNum * m_camsNumInStrip;

  intrins.resize(camNum, Intrin(m_f / m_pixSize));
  extrins.resize(camNum);
  images.resize(camNum);
  pts.resize(m_ptsNum);

  Scalar sceneXDim = m_imgW * m_grdRes + 
       lateralRangePerCam * (m_stripsNum - 1);
  Scalar sceneYDim = m_imgH * m_grdRes +
       longitudinalRangePerCam * (m_camsNumInStrip - 1);
  Scalar sceneZDim = m_sceneMaxHeight;
  Mat33 camPosCov = Mat33::Identity();
  camPosCov(0, 0) = m_camPlannarDev;
  camPosCov(1, 1) = m_camPlannarDev;
  camPosCov(2, 2) = m_camHeightDev;
  Mat33 camRotCov = Mat33::Identity();
  camRotCov *= m_camRotDev;
  Scalar scale = lateralRangePerCam;
  Scalar nwAngleRad = m_nwAngle / 180 * Scalar(M_PI);
  Mat33 nwRot;
  nwRot << cos(nwAngleRad), -sin(nwAngleRad), 0,
     sin(nwAngleRad),  cos(nwAngleRad), 0,
     0, 0, 1;
  for (size_t i = 0; i < m_stripsNum; i++)
  {
    for (size_t j = 0; j < m_camsNumInStrip; j++)
    {
      size_t id = i * m_camsNumInStrip + j;
      images[id].m_width = m_imgW;
      images[id].m_height = m_imgH;
      images[id].m_id = id;
      images[id].m_path = "test";

      Pos meanCamPos;
      meanCamPos << -sceneXDim * 0.5 + m_imgW * m_grdRes * 0.5 +
          i * lateralRangePerCam,
          -sceneYDim * 0.5 + m_imgH * m_grdRes * 0.5 + 
          j * longitudinalRangePerCam,
          flightHeight;

      hs::math::random::NormalRandomVar<Scalar, 3>::normRandomVar(
      meanCamPos, camPosCov, extrins[id].m_c);
      extrins[id].m_c = nwRot * extrins[id].m_c;
      //extrins[id].m_c /= scale;
      Pos meanAngles = Pos::Zero();
      Pos angles;
      hs::math::random::NormalRandomVar<Scalar, 3>::normRandomVar(
      meanAngles, camRotCov, angles);

      hs::math::geometry::EulerAngles<Scalar> ea(angles[0] / 180 * Scalar(M_PI),
             angles[1] / 180 * Scalar(M_PI),
             angles[2] / 180 * Scalar(M_PI));
      Mat33 R = ea.template toOrthoRotMat<2, 1, -3, 1>();
      R.transposeInPlace();
      extrins[id].m_r = R * nwRot.transpose();
    }
  }

  Pt3D maxPt;
  maxPt << sceneXDim * 0.5,
     sceneYDim * 0.5,
     sceneZDim * 0.5;
  Pt3D minPt = -maxPt;

  for (size_t i = 0; i < m_ptsNum; i++)
  {
  hs::math::random::UniformRandomVar<Scalar, 3>::uniformRandomVar(
    minPt, maxPt, pts[i]);
  pts[i] = nwRot * pts[i];
  //pts[i] /= scale;
  }

  return 0;
}

private:
/**
 *  相机的焦距米为单位
 */
Scalar m_f;
/**
 *  航带数
 */
size_t m_stripsNum;
/**
 *  每条航带包含的相机数
 */
size_t m_camsNumInStrip;
/**
 *  影像的地面分辨率
 */
Scalar m_grdRes;
/**
 *  影像的宽度，像素为单位
 */
ImgDim m_imgW;
/**
 *  影像的高度，像素为单位
 */
ImgDim m_imgH;
/**
 *  像素大小，米为单位
 */
Scalar m_pixSize;
/**
 *  场景中空间点的数量
 */
size_t m_ptsNum;
/**
 *  旁向重叠率
 */
Scalar m_lateralOverlap;
/**
 *  航向重叠率
 */
Scalar m_longitudinalOverlap;
/**
 *  场景中最高点的高度
 */
Scalar m_sceneMaxHeight;
/**
 *  相机高度分布标准差
 */
Scalar m_camHeightDev;
/**
 *  相机平面方向分布标准差
 */
Scalar m_camPlannarDev;
/**
 *  相机旋转角标准差，度为单位
 */
Scalar m_camRotDev;
/**
 *  飞行方向的北偏东角，度为单位
 */
Scalar m_nwAngle;
};

template <typename _Scalar, typename _ImgDim>
class KeysGenerator
{
public:
typedef _Scalar Scalar;
typedef _ImgDim ImgDim;

typedef hs::sfm::IntrinParam<Scalar> Intrin;
typedef EIGEN_VECTOR(Intrin) IntrinContainer;
typedef hs::sfm::ExtrinParam<Scalar> Extrin;
typedef typename Extrin::Pos Pos;
typedef typename Extrin::Rot Rot;
typedef EIGEN_VECTOR(Extrin) ExtrinContainer;
typedef EIGEN_VEC(Scalar, 3) Pt3D;
typedef EIGEN_VECTOR(Pt3D) Pt3DContainer;
typedef EIGEN_VEC(Scalar, 2) Pt2D;
typedef EIGEN_VECTOR(Pt2D) Pt2DSet;
typedef EIGEN_VECTOR(Pt2DSet) Pt2DSetContainer;
typedef EIGEN_MAT(Scalar, 3, 3) Mat33;
typedef CamFunc<Scalar> Cam;
typedef typename Cam::PMat PMat;
typedef ImageKeys<Scalar> ImgKeys;
typedef EIGEN_VECTOR(ImgKeys) ImgKeysContainer;
typedef std::vector<std::pair<size_t, size_t> > Track;
typedef std::vector<Track> TrackContainer;

typedef int Err;

KeysGenerator(ImgDim imgW, ImgDim imgH) : m_imgW(imgW), m_imgH(imgH) {}

Err operator()(const IntrinContainer& intrins, 
     const ExtrinContainer& extrins, 
     const Pt3DContainer& pts,
     ImgKeysContainer& imgKeysSet, 
     TrackContainer& tracks) const
{
  size_t camNum = intrins.size();
  if (camNum != extrins.size())
  {
  return -1;
  }

  Pt2DSetContainer keySets(camNum);
  size_t ptsNum = pts.size();
  tracks.resize(ptsNum);
  for (size_t i = 0; i < ptsNum; i++)
  {
    const Pt3D& pt = pts[i];
    for (size_t j = 0; j < camNum; j++)
    {
      PMat P = Cam::getPMat(intrins[j], extrins[j]);
      Pt3D keyH = P.block(0, 0, 3, 3) * pts[i] + 
        P.block(0, 3, 3, 1);
      keyH /= keyH(2);
      if (keyH(0) > (-Scalar(m_imgW) / 2) && 
      keyH(0) < ( Scalar(m_imgW) / 2) &&
      keyH(1) > (-Scalar(m_imgH) / 2) && 
      keyH(1) < ( Scalar(m_imgH) / 2))
      {
      tracks[i].push_back(
        std::make_pair(j, keySets[j].size()));
      keySets[j].push_back(keyH.segment(0, 2));
      }
    }
  }
  imgKeysSet.clear();
  for (size_t i = 0; i < camNum; i++)
  {
  imgKeysSet.push_back(ImgKeys(keySets[i]));
  }

  return 0;
}

private:
ImgDim m_imgW;
ImgDim m_imgH;
};

class MatchGenerator
{
public:
typedef std::vector<std::pair<size_t, size_t> > Track;
typedef std::vector<Track> TrackContainer;
typedef int Err;

Err operator()(const TrackContainer& tracks,
     hs::sfm::MatchContainer& matches) const
{
  using namespace hs::sfm;
  size_t tracksNum = tracks.size();
  matches.clear();
  for (size_t i = 0; i < tracksNum; i++)
  {
  size_t viewNum = tracks[i].size();
  for (size_t j = 0; j < viewNum; j++)
  {
    for (size_t k = j + 1; k < viewNum; k++)
    {
    ImgPair imgPair(tracks[i][j].first, tracks[i][k].first);
    
    matches[imgPair].push_back(
      KeyPair(tracks[i][j].second, tracks[i][k].second));
    }
  }
  }

  return 0;
}
};

template <typename _Scalar>
struct NoiseImageKeys
{
typedef _Scalar Scalar;
typedef ImageKeys<Scalar> ImgKeys;
typedef EIGEN_VEC(Scalar, 2) Key;
typedef EIGEN_MAT(Scalar, 2, 2) Cov;
typedef int Err;

NoiseImageKeys(const Cov& cov)
  : m_cov(cov) {}

Err operator()(ImgKeys& imgKeys) const
{
  size_t imgKeysNum = imgKeys.size();
  for (size_t i = 0; i < imgKeysNum; i++)
  {
  Key mean = imgKeys[i];
  hs::math::random::NormalRandomVar<Scalar, 2>::normRandomVar(
    mean, m_cov, imgKeys[i]);
  }

  return 0;
}

Cov m_cov;
};

template <typename _Scalar>
struct NoiseExtrin
{
typedef _Scalar Scalar;
typedef ExtrinParam<Scalar> Extrin;
typedef typename Extrin::Pos Pos;
typedef typename Extrin::Rot Rot;
typedef hs::math::geometry::EulerAngles<Scalar> EA;
typedef typename EA::OrthoRotMat RMat;
typedef EIGEN_MAT(Scalar, 3, 3) Cov;
typedef int Err;

NoiseExtrin(const Cov& covR, const Cov& covC)
  : m_covR(covR), m_covC(covC) {}

Err operator()(Extrin& extrin) const
{
  Pos meanC = extrin.m_c;
  hs::math::random::NormalRandomVar<Scalar, 3>::normRandomVar(
  meanC, m_covC, extrin.m_c);

  RMat R = extrin.m_r;
  R.transposeInPlace();
  EA ea;
  ea.template fromOrthoRotMat<2, 1, -3, 1>(R);
  Pos meanR;
  meanR << ea[0],
  ea[1],
  ea[2];
  Pos noiseAngles;
  hs::math::random::NormalRandomVar<Scalar, 3>::normRandomVar(
  meanR, m_covR, noiseAngles);
  ea[0] = noiseAngles[0];
  ea[1] = noiseAngles[1];
  ea[2] = noiseAngles[2];
  R = ea.template toOrthoRotMat<2, 1, -3, 1>();
  R.transposeInPlace();
  extrin.m_r = R;

  return 0;
}

Cov m_covR;
Cov m_covC;
};

}//sfm
}//hs

#endif
