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
 *  ��򵥵�BA��������������������ά����ΪX����������ΪY��
 *  ��ά�����������Ϊ��������꣬����������ڲ�����
 *  �����תʹ�������ת��
 */
template <typename _Scalar>
class BANaiveVecFunc
{
public:
  typedef _Scalar Scalar;

  /**
   *  X�������������������ά�����꣬�����������ǰ��
   */
  typedef EIGEN_VEC(Scalar, Eigen::Dynamic) XVec;
  /**
   *  Y����������ά����������
   */
  typedef EIGEN_VEC(Scalar, Eigen::Dynamic) YVec;

  typedef EIGEN_VEC(Scalar, 3) Vec3;

  typedef typename XVec::Index Index;

  typedef int Err;

  /**
   *  ÿ̨����Ĳ�������Ϊ6��������תռ3����λ��ռ3��
   */
  static const Index m_paramsPerCam = 6;

  /**
   *  ��ά���������Ϊ3
   */
  static const Index m_paramsPerPt = 3;

  /**
   *  �������������Ϊ2
   */
  static const Index m_paramsPerFeat = 2;

  /**
   *  ��ʾ���������������ά��֮���ӳ���ϵ
   *  pair��firstֵ��ʾ�����ţ�secondֵ��ʾ��ά�����
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

      Vec3 pt = x.segment(camParamSize + ptId * m_paramsPerPt,
                m_paramsPerPt);

      Vec3 r = x.segment(camId * m_paramsPerCam, 3);
      Vec3 t = x.segment(camId * m_paramsPerCam + 3, 3);

      //Rot3D<Scalar> R(r);
      //Vec3 ptCam = R * pt + t;

      //�����ת����ƽ��
      Scalar theta = r.norm();
      r /= theta;
      Vec3 ptCam = cos(theta) * pt + 
        sin(theta) * r.cross(pt) + 
        (1 - cos(theta)) * pt.dot(r) * r + t;

      ptCam /= ptCam[2];

      y.segment(i * m_paramsPerFeat, m_paramsPerFeat) = 
        ptCam.segment(0, 2);
    }

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
