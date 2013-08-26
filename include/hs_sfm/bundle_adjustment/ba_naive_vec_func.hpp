#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_VEC_FUNC_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_VEC_FUNC_HPP_

#include <vector>
#include <utility>

#include "hs_math/linear_algebra/eigen_macro.hpp"
//#include "hs_math/geometry/rotation.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

/**
 *  最简单的BA函数。将相机外参数与三维点作为X，特征点作为Y。
 *  三维点与特征点均为非齐次坐标，不考虑相机内参数。
 *  相机旋转使用轴角旋转。
 */
template <typename _Scalar>
class BANaiveVecFunc
{
public:
  typedef _Scalar Scalar;

  /**
   *  X向量包含相机参数与三维点坐标，相机参数排在前边
   */
  typedef EIGEN_VEC(Scalar, Eigen::Dynamic) XVec;
  /**
   *  Y向量包含二维特征点坐标
   */
  typedef EIGEN_VEC(Scalar, Eigen::Dynamic) YVec;

  typedef EIGEN_VEC(Scalar, 3) Vec3;
  typedef EIGEN_VEC(Scalar, 2) Vec2;

  typedef typename XVec::Index Index;

  typedef int Err;

  /**
   *  每台相机的参数个数为6，其中旋转占3个，位置占3个
   */
  static const Index m_paramsPerCam = 6;

  /**
   *  三维点参数个数为3
   */
  static const Index m_paramsPerPt = 3;

  /**
   *  特征点参数个数为2
   */
  static const Index m_paramsPerFeat = 2;

  /**
   *  表示特征点与相机、三维点之间的映射关系
   *  pair的first值表示相机序号，second值表示三维点序号
   */
  typedef std::pair<Index, Index> FeatMap;
  typedef std::vector<FeatMap> FeatMapContainer;

public:
  BANaiveVecFunc() {}
  BANaiveVecFunc(Index camNum, Index ptNum, Index featNum,
           const FeatMapContainer& featMaps)
    : m_camNum(camNum), m_ptNum(ptNum), m_featNum(featNum),
      m_featMaps(featMaps) {}

  inline Index getCamNum() const
  {
    return m_camNum;
  }

  inline void setCamNum(Index camNum)
  {
    m_camNum = camNum;
  }

  inline Index getPtNum() const
  {
    return m_ptNum;
  }

  inline void setPtNum(Index ptNum)
  {
    m_ptNum = ptNum;
  }

  inline Index getFeatNum() const
  {
    return m_featNum;
  }

  inline void setFeatNum(Index featNum)
  {
    m_featNum = featNum;
  }

  inline const FeatMapContainer& getFeatMaps() const
  {
    return m_featMaps;
  }

  inline void setFeatMaps(const FeatMapContainer& featMaps)
  {
    m_featMaps = featMaps;
  }

  inline Index getXSize() const
  {
    return getCamParamsSize() + getPtParamsSize();
  }

  inline Index getYSize() const
  {
    return m_featNum * m_paramsPerFeat;
  }

  inline Index getCamParamsSize() const
  {
    return m_camNum * m_paramsPerCam;
  }

  inline Index getPtParamsSize() const
  {
    return m_ptNum * m_paramsPerPt;
  }

  Err operator() (const XVec& x, YVec& y) const
  {
    Index xSize = x.rows();
    if (xSize != getXSize()) return -1;
    if (m_featNum != Index(m_featMaps.size())) return -1;
    Index ySize = getYSize();
    y.resize(ySize);

    Index camParamSize = getCamParamsSize();

    for (Index i = 0; i < m_featNum; i++)
    {
      Index camId = m_featMaps[i].first;
      Index ptId = m_featMaps[i].second;

      Vec3 p = x.segment(camParamSize + ptId * m_paramsPerPt,
               m_paramsPerPt);

      Vec3 r = x.segment(camId * m_paramsPerCam, 3);
      Vec3 t = x.segment(camId * m_paramsPerCam + 3, 3);

      Vec2 feat;
      ptPrjToFeat(r, t, p, feat);

      y.segment(i * m_paramsPerFeat, m_paramsPerFeat) = feat;
    }

    return 0;
  }

  inline static Err ptPrjToFeat(const Vec3& r, const Vec3& t, const Vec3& p,
                                Vec2& feat)
  {
    //hs::math::geometry::Rot3D<Scalar> R(r);
    //Vec3 ptCam = R * p + t;

    //轴角旋转加上平移
    Scalar theta = r.norm();
    Vec3 rN = r / theta;
    Vec3 ptCam = cos(theta) * p +
      sin(theta) * rN.cross(p) +
      (1 - cos(theta)) * p.dot(rN) * rN + t;
    ptCam /= ptCam[2];
    feat = ptCam.segment(0, 2);

    return 0;
  }

private:
  Index m_camNum;
  Index m_ptNum;
  Index m_featNum;
  FeatMapContainer m_featMaps;
};

}//ba
}//sfm
}//hs

#endif
