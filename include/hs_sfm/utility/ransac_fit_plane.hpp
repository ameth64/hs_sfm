#ifndef _HS_SFM_UTILITY_RANSAC_FIT_PLANE_HPP_
#define _HS_SFM_UTILITY_RANSAC_FIT_PLANE_HPP_

#include "hs_fit/ransac/ransac.hpp"
#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{

template <typename _Pt>
struct PtScalarType
{
  typedef typename _Pt::Scalar type;
};

//template <typename Scalar>
//struct PtScalarType<Scalar[3]>
//{
//  typedef Scalar type;
//};

//template <typename Scalar>
//struct PtScalarType<Scalar*>
//{
//  typedef Scalar type;
//};

template <typename _Pt>
class RansacFitPlane
{
public:
  typedef _Pt Pt;
  typedef typename hs::fit::PtSetType<Pt>::type PtSet;
  typedef typename PtScalarType<Pt>::type Scalar;
  EIGEN_BASE_DEF(Scalar);
  typedef Vec4 Plane;

  typedef int Err;

  //计算给定点集的正交回归
  class PlaneOrthoRegCalc
  {
  public:
    typedef Plane RansacModel;
    static const size_t m_modelMinSetSize = 3;

    int operator()(const PtSet& pts, Plane& pln)
    {
      size_t ptsNum = pts.size();
      if (ptsNum < m_modelMinSetSize)
      {
        return -1;
      }

      //pca
      //求均值
      Vec3 mean = Vec3::Zero();
      for (size_t i = 0; i < ptsNum; i++)
      {
        mean[0] += pts[i][0];
        mean[1] += pts[i][1];
        mean[2] += pts[i][2];
      }
      mean /= Scalar(ptsNum);
      Mat3X A(3, ptsNum);
      for (size_t i = 0; i < ptsNum; i++)
      {
        Vec3 pt;
        pt << pts[i][0],
          pts[i][1],
          pts[i][2];
        A.col(i) = pt - mean;
      }
      //协方差矩阵
      Mat33 cov = A * A.transpose();
      //求协方差矩阵最小特征值对应的特征向量
      Eigen::EigenSolver<Mat33> es(cov);
      Scalar minEigenVal = std::numeric_limits<Scalar>::max();
      int minIdx = -1;
      for (int i = 0; i < 3; i++)
      {
        Scalar eigenVal = es.eigenvalues()[i].real();
        if (minEigenVal > eigenVal)
        {
          minEigenVal = eigenVal;
          minIdx = i;
        }
      }
      Vec3 plnNorm = es.eigenvectors().col(minIdx).real();
      plnNorm /= plnNorm.norm();

      //计算过中值点，且法向量为pca方向的平面
      pln.segment(0, 3) = plnNorm;
      pln[3] = -plnNorm.dot(mean);

      return 0;
    }
  };

  class PtPlaneDistCalc
  {
  public:
    typedef Scalar Distance;

    int operator()(const Pt& pt, const Plane& pln, Distance& dist)
    {
      Distance norm = pln.segment(0, 3).norm();
      if (norm == 0)
      {
        return -1;
      }

      dist = std::abs(pt[0] * pln[0] + 
              pt[1] * pln[1] + 
              pt[2] * pln[2] + pln[3]) / norm;

      return 0;
    }
  };

  typedef hs::fit::Ransac<Pt, PlaneOrthoRegCalc, PtPlaneDistCalc> FitRansac;
  typedef typename FitRansac::IndexSet IndexSet;

  Err operator()(const PtSet& pts, Plane& pln, Scalar distThd) const
  {
    PlaneOrthoRegCalc modelCalc;
    PtPlaneDistCalc distCalc;
    FitRansac fitRansac(modelCalc, distCalc);

    fitRansac.setAlphaThd(0.95);
    fitRansac.setDistThd(distThd);
    PtSet refinedPts;
    IndexSet inlierIdx;
    PtSet bestPts;
    if (fitRansac(pts, refinedPts, inlierIdx, bestPts) != 0)
    {
      return -1;
    }

    return modelCalc(refinedPts, pln);
  }
};

}//sfm
}//hs

#endif
