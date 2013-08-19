#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_ANALYTIC_JAC_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_ANALYTIC_JAC_HPP_

#include <vector>

#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_jac_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VecFunc>
class BANaiveAnalyticJacobian;

/**
 *  计算简单Bundle Adjustment函数的Jacbobian矩阵。
 *
 *  设\f$\mathbf{y_i}\f$为第i个影像特征点对应\f$\mathbf{c_j}\f$为第j个相机的参数，
 *  \f$\mathbf{p_k}\f$为第k个三维点，
 *  则\f$\mathbf{y_i}\f$的齐次表示\f$\mathbf{\hat{y}_i}\f$有：
 *  \f[
 *    \mathbf{\hat{y}_i} = \mathbf{r_j}(\mathbf{p_k}) + \mathbf{t_j}
 *  \f]
 *  其中\f$\mathbf{t_j}\f$表示第j个相机的位移。
 *  \f$\mathbf{r_j}(\mathbf{x})\f$表示第j个相机的朝向，使用轴角旋转表示为：
 *  \f[
 *    \mathbf{r_j}(\mathbf{x}) = cos(\theta)\mathbf{x} +
 *      sin(\theta)((\mathbf{r_j} / \theta) \times \mathbf{x}) +
 *      \frac{1-cos(\theta)}{\theta^2}(\mathbf{r_j}\cdot\mathbf{x})\mathbf{r_j}
 *  \f]
 *  其中\f$\theta\f$为\f$\mathbf{r_j}\f$的模。
 *  于是有\f$\mathbf{\hat{y}_i} = [\hat{y}_{i0},\hat{y}_{i1},\hat{y}_{i2}]^T\f$：
 *  \f[
 *    \hat{y}_{i0} = cos(\theta)p_{k0} + 
 *           \frac{sin(\theta)}{\theta}r_{j1}p_{k2} +
 *           d\frac{1-cos(\theta)}{\theta^2}r_{j0} + t_{j0}
 *  \f]
 *  \f[
 *    \hat{y}_{i1} = cos(\theta)p_{k1} + 
 *           \frac{sin(\theta)}{\theta}r_{j2}p_{k0} +
 *           d\frac{1-cos(\theta)}{\theta^2}r_{j1} + t_{j1}
 *  \f]
 *  \f[
 *    \hat{y}_{i2} = cos(\theta)p_{k2} + 
 *           \frac{sin(\theta)}{\theta}r_{j0}p_{k1} +
 *           d\frac{1-cos(\theta)}{\theta^2}r_{j2} + t_{j2}
 *  \f]
 *  其中\f$d=\mathbf{r_j}\cdot\mathbf{p_k}\f$
 */
template <typename _Scalar>
class BANaiveAnalyticJacobian<BANaiveVecFunc<_Scalar> >
{
public:
  typedef _Scalar Scalar;

  typedef int Err;

  typedef BANaiveVecFunc<Scalar> VecFunc;
  typedef typename VecFunc::XVec XVec;
  typedef typename VecFunc::YVec YVec;
  typedef typename VecFunc::Index Index;
  typedef typename VecFunc::FeatMap FeatMap;
  typedef typename VecFunc::FeatMapContainer FeatMapContainer;
  typedef BANaiveJacMatrix<Scalar, Index,
                           VecFunc::m_paramsPerFeat,
                           VecFunc::m_paramsPerCam,
                           VecFunc::m_paramsPerPt> Jac;
  typedef typename Jac::DrvIdx DrvIdx;

  Err operator()(const VecFunc& vecFunc, const XVec& x, Jac& jac) const
  {
    typedef EIGEN_VEC(Scalar, 3) Vec3;
    typedef EIGEN_MAT(Scalar, 3, 3) Mat33;
    typedef Eigen::Triplet<DrvIdx, Index> TripletType;
    Index camNum = vecFunc.getCamNum();
    Index ptNum = vecFunc.getPtNum();
    Index featNum = vecFunc.getFeatNum();
    Index camParamsSize = vecFunc.getCamParamsSize();
    FeatMapContainer featMaps = vecFunc.getFeatMaps();
    Index xSize = vecFunc.getXSize();
    std::vector<TripletType> camsTriplet;
    std::vector<TripletType> ptsTriplet;
    jac.clear();
    jac.m_camNum = camNum;
    jac.m_ptNum = ptNum;
    for (Index i = 0; i < featNum; i++)
    {
      const FeatMap& featMap = featMaps[i];
      Index j = featMap.first;
      Index k = featMap.second;

      Vec3 r = x.segment(j * VecFunc::m_paramsPerCam, 3);
      Vec3 t = x.segment(j * VecFunc::m_paramsPerCam + 3, 3);
      Vec3 p = x.segment(camParamsSize + k * VecFunc::m_paramsPerPt, 
                 VecFunc::m_paramsPerPt);

      Scalar theta = r.norm();
      Scalar theta2 = theta * theta;
      Scalar theta4 = theta2 * theta2;
      Scalar d = r.dot(p);

      //sin(theta)
      Scalar st = sin(theta);
      //cos(theta)
      Scalar ct = cos(theta);

      //f = cos(theta) * p
      Vec3 f = ct * p;
      //\frac{\partial f}{\partial r}
      Mat33 pfpr = p * (-st / theta * r).transpose();
      //\frac{\partial f}{\partial p}
      Mat33 pfpp = Mat33::Identity() * ct;

      //s = sin(theta) / theta
      Scalar s = st / theta;
      //\frac{\partial s}{\partial r}
      Vec3 pspr = (ct * theta - st) / theta2 / theta * r;

      //v = r x p
      Vec3 v = r.cross(p);
      //\frac{\partial v}{\partial r}
      Mat33 pvpr;
      pvpr << 0, p[2], -p[1],
          -p[2], 0, p[0],
          p[1], -p[0], 0;
      //\frac{\partial v}{\partial p}
      Mat33 pvpp;
      pvpp << 0, -r[2], r[1],
          r[2], 0, -r[0],
          -r[1], r[0], 0;

      //c = (1 - cos(theta)) / theta^2
      Scalar c = (1 - ct) / theta2;
      //\frac{partial c}{\partial r}
      Vec3 pcpr = 
        (st * theta - 2 * (1 - ct)) / theta4 * r;

      //u = r^T * p * r
      Vec3 u = d * r;
      //\frac{partial u}{\partial r}
      Mat33 pupr;
      pupr << d + r[0] * p[0], r[0] * p[1], r[0] * p[2],
          r[1] * p[0], d + r[1] * p[1], r[1] * p[2],
          r[2] * p[0], r[2] * p[1], d + r[2] * p[2];
      //\frac{partial u}{\patial p}
      Mat33 pupp;
      pupp << r[0] * r[0], r[0] * r[1], r[0] * r[2],
          r[1] * r[0], r[1] * r[1], r[1] * r[2],
          r[2] * r[0], r[2] * r[1], r[2] * r[2];

      Vec3 y = f + s * v + c * u + t;
      //\frac{partial y}{\partial r}
      Mat33 pypr = pfpr +
             v * pspr.transpose() + s * pvpr +
             u * pcpr.transpose() + c * pupr;
      //\frac{partial y}{\partial p}
      Mat33 pypp = pfpp + s * pvpp + c * pupp;
      //\frac{partial y}{\partial t} is identity
      Mat33 pypt = Mat33::Identity();

      typename Jac::CamDrvBlk camDrvBlk;
      camDrvBlk.m_camId = j;
      camDrvBlk.m_featId = i;
      for (Index m = 0; m < 3; m++)
      {
        //for r
        camDrvBlk.m_drvBlk(0, m) =
          (pypr(0, m) * y[2] - y[0] * pypr(2, m)) / y[2] / y[2];
        camDrvBlk.m_drvBlk(1, m) =
          (pypr(1, m) * y[2] - y[1] * pypr(2, m)) / y[2] / y[2];

        //for t
        camDrvBlk.m_drvBlk(0, 3 + m) =
          (pypt(0, m) * y[2] - y[0] * pypt(2, m)) / y[2] / y[2];
        camDrvBlk.m_drvBlk(1, 3 + m) =
          (pypt(1, m) * y[2] - y[1] * pypt(2, m)) / y[2] / y[2];
      }

      jac.m_camsDrv.push_back(camDrvBlk);
      camsTriplet.push_back(TripletType(j, k, jac.m_camsDrv.size()));

      typename Jac::PtDrvBlk ptDrvBlk;
      ptDrvBlk.m_ptId = k;
      ptDrvBlk.m_featId = i;
      for (Index m = 0; m < 3; m++)
      {
        ptDrvBlk.m_drvBlk(0, m) =
          (pypp(0, m) * y[2] - y[0] * pypp(2, m)) / y[2] / y[2];
        ptDrvBlk.m_drvBlk(1, m) =
          (pypp(1, m) * y[2] - y[1] * pypp(2, m)) / y[2] / y[2];
      }

      jac.m_ptsDrv.push_back(ptDrvBlk);
      ptsTriplet.push_back(TripletType(j, k, jac.m_ptsDrv.size()));
    }

    jac.m_camsDrvMap.resize(camNum, ptNum);
    jac.m_camsDrvMap.setFromTriplets(camsTriplet.begin(), 
                     camsTriplet.end());
    jac.m_ptsDrvMap.resize(camNum, ptNum);
    jac.m_ptsDrvMap.setFromTriplets(ptsTriplet.begin(),
                    ptsTriplet.end());

    return 0;
  }

};

}//ba
}//sfm
}//hs

#endif
