#ifndef _HS_SFM_UTILITY_SFM_FILE_IO_HPP_
#define _HS_SFM_UTILITY_SFM_FILE_IO_HPP_

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/ransac_fit_plane.hpp"

namespace hs
{
namespace fit
{

template <typename _Scalar>
struct PointSetType<EIGEN_VECTOR(_Scalar, 3) >
{
  typedef EIGEN_VECTOR(_Scalar, 3) Pt;
  typedef EIGEN_STD_VECTOR(Pt) type;
};

}

namespace sfm
{

template <typename _Scalar>
struct LoadImageKeys
{
  typedef _Scalar Scalar;
  typedef ImageKeys<Scalar> Keys;
  typedef int Err;

  Err operator()(const std::string& keys_path,
                 Keys& keys) const
  {
    std::ifstream keys_file(keys_path.c_str(), std::ios::in);
    if (!keys_file.is_open())
    {
      return -1;
    }

    size_t number_of_keys;
    keys_file>>number_of_keys;
    keys.resize(number_of_keys);
    for (size_t i = 0; i < number_of_keys; i++)
    {
      keys_file>>keys[i][0]>>keys[i][1];
    }

    return 0;
  }
};

template <typename _Scalar>
struct SaveImageKeys
{
  typedef _Scalar Scalar;
  typedef ImageKeys<Scalar> Keys;
  typedef int Err;
  
  Err operator()(const std::string& keys_path, 
                 const Keys& keys) const
  {
    std::ofstream keys_file(keys_path.c_str(), std::ios::out);
    if (!keys_file.is_open())
    {
      return -1;
    }

    keys_file.setf(std::ios::fixed);
    keys_file<<std::setprecision(6);

    size_t number_of_keys = keys.size();
    keys_file<<number_of_keys<<'\n';
    for (size_t i = 0; i < number_of_keys; i++)
    {
      keys_file<<keys[i][0]<<' '<<keys[i][1]<<'\n';
    }

    return 0;
  }
};

struct LoadMatches
{
  typedef hs::sfm::ImagePair ImagePair;
  typedef hs::sfm::KeyPair KeyPair;
  typedef hs::sfm::KeyPairContainer KeyPairContainer;
  typedef hs::sfm::MatchContainer MatchContainer;
  typedef int Err;
  Err operator()(const std::string& matches_path,
                 MatchContainer& matches) const
  {
    std::ifstream matches_file(matches_path.c_str(), std::ios::in);
    if (!matches_file.is_open())
    {
      return -1;
    }

    matches.clear();
    while (!matches_file.eof())
    {
      char buf[512];
      matches_file.getline(buf, 512);
      if (buf[0] == 0)
      {
        break;
      }
      std::stringstream ss(buf);
      ImagePair image_pair;
      ss>>image_pair.first>>image_pair.second;
      size_t number_of_image_pairs;
      matches_file.getline(buf, 512);
      ss.clear();
      ss.str(buf);
      ss>>number_of_image_pairs;
      KeyPairContainer key_pairs(number_of_image_pairs);
      for (size_t i = 0; i < number_of_image_pairs; i++)
      {
        matches_file.getline(buf, 512);
        ss.clear();
        ss.str(buf);
        ss>>key_pairs[i].first>>key_pairs[i].second;
      }
      matches[image_pair] = key_pairs;
    }

    return 0;
  }
};

struct SaveMatches
{
  typedef hs::sfm::MatchContainer Matches;
  typedef int Err;
  Err operator()(const std::string& matches_path,
                 const Matches& matches) const
  {
    std::ofstream matches_file(matches_path.c_str(), std::ios::out);
    if (!matches_file.is_open())
    {
      return -1;
    }

    Matches::const_iterator itr = matches.begin();
    Matches::const_iterator itr_end = matches.end();
    for (; itr != itr_end; ++itr)
    {
      matches_file<<(itr->first).first<<" "<<(itr->first).second<<"\n";
      size_t number_of_key_matches = itr->second.size();
      matches_file<<number_of_key_matches<<"\n";
      for (size_t i = 0; i < number_of_key_matches; i++)
      {
        matches_file<<itr->second[i].first<<" "<<
                      itr->second[i].second<<"\n";
      }
    }

    matches_file.close();

    return 0;
  }
};

template <typename _Scalar, typename _ImageDimension>
struct SaveXugFile
{
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;

  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinContainer;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef typename ExtrinsicParams::Position Position;
  typedef typename ExtrinsicParams::Rotation Rotation;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinContainer;
  typedef ImageParams<ImageDimension> Image;
  typedef EIGEN_STD_VECTOR(Image) ImageContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;
  typedef ImageKeys<Scalar> Keys;
  typedef EIGEN_STD_VECTOR(Keys) KeysContainer;
  typedef std::vector<std::pair<size_t, size_t> > Track;
  typedef std::vector<Track> TrackContainer;
  typedef RansacFitPlane<Point3D> PlaneFitter;
  typedef typename PlaneFitter::Plane Plane;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

  typedef int Err;
  
  Err operator() (const std::string& xug_path,
                  const IntrinContainer& intrinsic_params_set, 
                  const ExtrinContainer& extrinsic_params_set,
                  const ImageContainer& images,
                  const Point3DContainer& points,
                  Scalar fit_threshold = Scalar(0.05)) const
  {
    size_t number_of_cameras = intrinsic_params_set.size();
    if (number_of_cameras != extrinsic_params_set.size() ||
        number_of_cameras != images.size())
    {
      return -1;
    }
    size_t number_of_points = points.size();

    std::ofstream xug_file(xug_path.c_str(), std::ios::out);
    if (!xug_file.is_open())
    {
      return -1;
    }
    xug_file<<"# Bundle file v0.3\n"; 
    xug_file<<number_of_cameras<<' '<<number_of_points<<'\n';

    //通过点云拟合场景平面
    PlaneFitter fitter;
    Plane plane;
    if (fitter(points, plane, fit_threshold) != 0) return -1;
    xug_file<<plane[0]<<' '<<plane[1]<<' '<<plane[2]<<' '<<plane[3]<<'\n';
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      xug_file<<images[i].m_path<<'\n';
      xug_file<<"0 0 0 0\n";
      Position c = extrinsic_params_set[i].position();
      Matrix33 R = extrinsic_params_set[i].rotation();
      //计算相机方向的四个顶点
      Scalar distance = std::abs(c.dot(plane.segment(0, 3)) +
                                 plane[3]) / plane.segment(0, 3).norm();
      Scalar depth = distance * 0.1;
      Scalar w = Scalar(images[i].m_width) / 
             intrinsic_params_set[i].focal_length() * depth * Scalar(0.5);
      Scalar h = Scalar(images[i].m_height) /
             intrinsic_params_set[i].focal_length() * depth * Scalar(0.5);
      Position corner[4];
      corner[0] << -w, -h, -depth;
      corner[1] <<  w, -h, -depth;
      corner[2] <<  w,  h, -depth;
      corner[3] << -w,  h, -depth;
      for (int j = 0; j < 4; j++)
      {
        corner[j] = R.inverse() * corner[j] + c;
        xug_file<<corner[j][0]<<' '<<corner[j][1]<<' '<<corner[j][2]
             <<'\n';
      }
      xug_file<<images[i].m_width<<' '<<images[i].m_height<<'\n';
      xug_file<<intrinsic_params_set[i].focal_length()<<" 0 0\n";
      xug_file<<R(0, 0)<<" "<<R(0, 1)<<" "<<R(0, 2)<<'\n'
              <<R(1, 0)<<" "<<R(1, 1)<<" "<<R(1, 2)<<'\n'
              <<R(2, 0)<<" "<<R(2, 1)<<" "<<R(2, 2)<<'\n';
      xug_file<<c[0]<<" "<<c[1]<<" "<<c[2]<<'\n';
    }

    //输出点云
    for (size_t i = 0; i < number_of_points; i++)
    {
      xug_file<<points[i][0]<<' '<<points[i][1]<<' '<<points[i][2]<<'\n';
      xug_file<<"255 255 0\n";
      xug_file<<"0\n";
    }

    return 0;
  }
};

}
}

#endif
