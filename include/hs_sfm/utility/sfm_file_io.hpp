#ifndef _HS_SFM_UTILITY_SFM_FILE_IO_HPP_
#define _HS_SFM_UTILITY_SFM_FILE_IO_HPP_

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "hs_sfm/utility/key_type.hpp"
#include "hs_sfm/utility/match_type.hpp"
#include "hs_sfm/utility/cam_type.hpp"
#include "hs_sfm/utility/ransac_fit_plane.hpp"

namespace hs
{
namespace fit
{

template <typename _Scalar>
struct PtSetType<EIGEN_VEC(_Scalar, 3) >
{
  typedef EIGEN_VEC(_Scalar, 3) Pt;
  typedef EIGEN_VECTOR(Pt) type;
};

}

namespace sfm
{

template <typename _Scalar>
struct LoadImageKeys
{
  typedef _Scalar Scalar;
  typedef ImageKeys<Scalar> ImgKeys;
  typedef int Err;

  Err operator()(const std::string& keysPath,
           ImgKeys& imgKeys) const
  {
    std::ifstream keyFile(keysPath.c_str(), std::ios::in);
    if (!keyFile.is_open())
    {
      return -1;
    }

    size_t keysNum;
    keyFile>>keysNum;
    imgKeys.resize(keysNum);
    for (size_t i = 0; i < keysNum; i++)
    {
      keyFile>>imgKeys[i][0]>>imgKeys[i][1];
    }

    return 0;
  }
};

template <typename _Scalar>
struct SaveImageKeys
{
  typedef _Scalar Scalar;
  typedef ImageKeys<Scalar> ImgKeys;
  typedef int Err;
  
  Err operator()(const std::string& keysPath, 
           const ImgKeys& imgKeys) const
  {
    std::ofstream keyFile(keysPath.c_str(), std::ios::out);
    if (!keyFile.is_open())
    {
      return -1;
    }

    keyFile.setf(std::ios::fixed);
    keyFile<<std::setprecision(6);

    size_t imgKeysNum = imgKeys.size();
    keyFile<<imgKeysNum<<'\n';
    for (size_t i = 0; i < imgKeysNum; i++)
    {
      keyFile<<imgKeys[i][0]<<' '<<imgKeys[i][1]<<'\n';
    }

    return 0;
  }
};

struct LoadMatches
{
  typedef hs::sfm::ImgPair ImgPair;
  typedef hs::sfm::KeyPair KeyPair;
  typedef hs::sfm::KeyPairContainer KeyPairContainer;
  typedef hs::sfm::MatchContainer MatchContainer;
  typedef int Err;
  Err operator()(const std::string& matchesPath,
           MatchContainer& matches) const
  {
    std::ifstream matchesFile(matchesPath.c_str(), std::ios::in);
    if (!matchesFile.is_open())
    {
      return -1;
    }

    matches.clear();
    while (!matchesFile.eof())
    {
      char buf[512];
      matchesFile.getline(buf, 512);
      if (buf[0] == 0)
      {
        break;
      }
      std::stringstream ss(buf);
      ImgPair imgPair;
      ss>>imgPair.first>>imgPair.second;
      size_t keyPairNum;
      matchesFile.getline(buf, 512);
      ss.clear();
      ss.str(buf);
      ss>>keyPairNum;
      KeyPairContainer keyPairs(keyPairNum);
      for (size_t i = 0; i < keyPairNum; i++)
      {
        matchesFile.getline(buf, 512);
        ss.clear();
        ss.str(buf);
        ss>>keyPairs[i].first>>keyPairs[i].second;
      }
      matches[imgPair] = keyPairs;
    }

    return 0;
  }
};

struct SaveMatches
{
  typedef hs::sfm::MatchContainer Matches;
  typedef int Err;
  Err operator()(const std::string& matchesPath,
           const Matches& matches) const
  {
    std::ofstream matchesFile(matchesPath.c_str(), std::ios::out);
    if (!matchesFile.is_open())
    {
      return -1;
    }

    Matches::const_iterator itr = matches.begin();
    Matches::const_iterator itrEnd = matches.end();
    for (; itr != itrEnd; ++itr)
    {
      matchesFile<<(itr->first).first<<" "<<(itr->first).second<<"\n";
      size_t keyMatchesCnt = itr->second.size();
      matchesFile<<keyMatchesCnt<<"\n";
      for (size_t i = 0; i < keyMatchesCnt; i++)
      {
        matchesFile<<itr->second[i].first<<" "<<
          itr->second[i].second<<"\n";
      }
    }

    matchesFile.close();

    return 0;
  }
};

template <typename _Scalar, typename _ImgDim>
struct SaveXugFile
{
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
  typedef ImageKeys<Scalar> ImgKeys;
  typedef EIGEN_VECTOR(ImgKeys) ImgKeysContainer;
  typedef std::vector<std::pair<size_t, size_t> > Track;
  typedef std::vector<Track> TrackContainer;
  typedef RansacFitPlane<Pt3D> PlaneFitter;
  typedef typename PlaneFitter::Plane Plane;
  typedef EIGEN_MAT(Scalar, 3, 3) Mat33;

  typedef int Err;
  
  Err operator() (const std::string& xugPath,
          const IntrinContainer& intrins, 
          const ExtrinContainer& extrins,
          const ImageContainer& images,
          const Pt3DContainer& pts,
          Scalar fitThd = Scalar(0.05)) const
  {
    size_t camNum = intrins.size();
    if (camNum != extrins.size() || camNum != images.size())
    {
      return -1;
    }
    size_t ptsNum = pts.size();

    std::ofstream xugFile(xugPath.c_str(), std::ios::out);
    if (!xugFile.is_open())
    {
      return -1;
    }
    xugFile<<"# Bundle file v0.3\n"; 
    xugFile<<camNum<<' '<<ptsNum<<'\n';

    //通过点云拟合场景平面
    PlaneFitter fit;
    Plane pln;
    if (fit(pts, pln, fitThd) != 0) return -1;
    xugFile<<pln[0]<<' '<<pln[1]<<' '<<pln[2]<<' '<<pln[3]<<'\n';
    for (size_t i = 0; i < camNum; i++)
    {
      xugFile<<images[i].m_path<<'\n';
      xugFile<<"0 0 0 0\n";
      Pos c = extrins[i].m_c;
      Mat33 R = extrins[i].m_r;
      //计算相机方向的四个顶点
      Scalar dist = std::abs(c.dot(pln.segment(0, 3)) + pln[3]) /
              pln.segment(0, 3).norm();
      Scalar depth = dist * 0.1;
      Scalar w = Scalar(images[i].m_width) / 
             intrins[i].m_focalLength * depth * Scalar(0.5);
      Scalar h = Scalar(images[i].m_height) /
             intrins[i].m_focalLength * depth * Scalar(0.5);
      Pos corner[4];
      corner[0] << -w, -h, -depth;
      corner[1] <<  w, -h, -depth;
      corner[2] <<  w,  h, -depth;
      corner[3] << -w,  h, -depth;
      for (int j = 0; j < 4; j++)
      {
        corner[j] = R.inverse() * corner[j] + c;
        xugFile<<corner[j][0]<<' '<<corner[j][1]<<' '<<corner[j][2]
             <<'\n';
      }
      xugFile<<images[i].m_width<<' '<<images[i].m_height<<'\n';
      xugFile<<intrins[i].m_focalLength<<" 0 0\n";
      xugFile<<R(0, 0)<<" "<<R(0, 1)<<" "<<R(0, 2)<<'\n'
           <<R(1, 0)<<" "<<R(1, 1)<<" "<<R(1, 2)<<'\n'
           <<R(2, 0)<<" "<<R(2, 1)<<" "<<R(2, 2)<<'\n';
      xugFile<<c[0]<<" "<<c[1]<<" "<<c[2]<<'\n';
    }

    //输出点云
    for (size_t i = 0; i < ptsNum; i++)
    {
      xugFile<<pts[i][0]<<' '<<pts[i][1]<<' '<<pts[i][2]<<'\n';
      xugFile<<"255 255 0\n";
      xugFile<<"0\n";
    }

    return 0;
  }
};

}
}

#endif
