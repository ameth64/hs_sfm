#ifndef _HS_SFM_UTILITY_CAM_TYPE_HPP_
#define _HS_SFM_UTILITY_CAM_TYPE_HPP_

#include <bitset>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_math/geometry/rotation.hpp"

namespace hs
{
namespace sfm
{

template <typename _Scalar>
class IntrinParam
{
public:
  typedef _Scalar Scalar;
  EIGEN_BASE_DEF(Scalar);
  /**
   *  内参数矩阵，表示相机内部的固有参数
   */
  typedef Mat33 KMat;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  IntrinParam()
    : m_focalLength(0), 
      m_skew(0),
      m_xPrinOffset(0),
      m_yPrinOffset(0),
      m_pixelRatio(1) {}
  IntrinParam(Scalar focalLength, 
        Scalar skew = 0, 
        Scalar xPrinOffset = 0,
        Scalar yPrinOffset = 0,
        Scalar pixelRatio = 1)
    : m_focalLength(focalLength),
      m_skew(skew),
      m_xPrinOffset(xPrinOffset),
      m_yPrinOffset(yPrinOffset),
      m_pixelRatio(pixelRatio) {}
  IntrinParam(const KMat& K)
    : m_focalLength(K(0, 0)), 
      m_skew(K(0, 1)),
      m_xPrinOffset(K(0, 2)),
      m_yPrinOffset(K(1, 2)),
      m_pixelRatio(K(0, 0) / K(1, 1)) {}

  inline KMat getKMat() const
  {
    KMat K;
    K << m_focalLength, m_skew, m_xPrinOffset,
       0, m_focalLength / m_pixelRatio, m_yPrinOffset,
       0, 0, 1;

    return K;
  }

  inline operator KMat() const
  {
    return getKMat();
  }

//private:
  /**
   *  焦距
   */
  Scalar m_focalLength;
  /**
   *  扭曲参数
   */
  Scalar m_skew;
  /**
   *  主点偏移
   */
  Scalar m_xPrinOffset;
  Scalar m_yPrinOffset;
  /**
   *  x方向与y方向上像素量纲比值
   */
  Scalar m_pixelRatio;
};

template <typename _Scalar>
class ExtrinParam
{
public:
  typedef _Scalar Scalar;
  EIGEN_BASE_DEF(Scalar);
  typedef hs::math::geometry::Rot3D<Scalar> Rot;
  typedef Vec3 Pos;
  /**
   *  外参数矩阵,表示相机运动
   */
  typedef Mat34 MMat;

  ExtrinParam() {m_r.setIdentity();m_c.setZero();}
  ExtrinParam(const Rot& r, const Pos& c)
    : m_r(r), m_c(c) {}

  inline MMat getMMat() const
  {
    MMat m;
    m.block(0, 0, 3, 3) = Mat33(m_r);
    m.block(0, 3, 3, 1) = -Mat33(m_r) * m_c;
    return m;
  }

  inline operator MMat() const
  {
    return getMMat();
  }

//private:
  /**
   *  相机朝向
   */
  Rot m_r;
  /**
   *  相机位置
   */
  Pos m_c;
};

template <typename _Dim>
class ImageParam
{
public:
  typedef _Dim Dim;
  
  size_t m_id;
  Dim m_width;
  Dim m_height;

  std::string m_path;
};

template <typename _Scalar>
class CamFunc
{
public:
  typedef _Scalar Scalar;
  EIGEN_BASE_DEF(Scalar);

  typedef IntrinParam<Scalar> Intrin;
  typedef ExtrinParam<Scalar> Extrin;

  typedef Mat34 PMat;
  typedef Mat33 RMat;

  static PMat getPMat(const Intrin& intr, const Extrin& extr)
  {
    PMat P = extr.getMMat();
    return intr.getKMat() * P;
  }
};

}
}

#endif
