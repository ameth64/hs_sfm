#include <iostream>
#include <ctime>

#include <gtest/gtest.h>

#include "hs_math/linear_algebra/latraits/vec_eigen.hpp"
#include "hs_math/linear_algebra/lafunc/arithmetic_eigen.hpp"
#include "hs_math/fdjac/ffd_sparse_jac.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_analytic_jac.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_ffd_jac.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_synthetic_data_generator.hpp"

#include "hs_math/fdjac/ffd_jac.hpp"

namespace
{

template <typename _Scalar>
class TestBANaiveJac
{
public:
  typedef _Scalar Scalar;

  typedef int Err;

  typedef hs::sfm::ba::BANaiveVecFunc<Scalar> BAVecFunc;
  typedef typename BAVecFunc::Index Index;
  typedef typename BAVecFunc::XVec XVec;
  typedef typename BAVecFunc::YVec YVec;

  typedef hs::math::fdjac::FwdFiniteDiffSparseJacobian<BAVecFunc> FFDJacobian;
  typedef typename FFDJacobian::Jac FFDJac;
  typedef hs::sfm::ba::BANaiveAnalyticJacobian<BAVecFunc> BAAnalyticJacobian;
  typedef typename BAAnalyticJacobian::Jac BAAnaliticJac;
  typedef hs::sfm::ba::BANaiveFFDJacobian<BAVecFunc> BAFFDJacobian;
  typedef typename BAFFDJacobian::Jac BAFFDJac;

  static Err test(const BAVecFunc& baVecFunc, const XVec& x, bool withFFD)
  {
    Index xSize = baVecFunc.getXSize();
    Index ySize = baVecFunc.getYSize();

    FFDJacobian ffdJac(Scalar(1e-6), Scalar(1e-8), Scalar(1e-10));
    BAAnalyticJacobian baAnalyticJac;
    BAFFDJacobian baFFDJac(Scalar(1e-6), Scalar(1e-8));

    FFDJac ffdJ;
    BAAnaliticJac baAnalyticJ;
    BAFFDJac baFFDJ;

    if (withFFD)
    {
      time_t ffdBegin = clock();
      if (0 != ffdJac(baVecFunc, x, ffdJ)) return -1;
      time_t ffdEnd = clock();
      std::cout<<"FFD took "
               <<Scalar(ffdEnd - ffdBegin) /
                 Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";
    }
    time_t baAnalyticBegin = clock();
    if (0 != baAnalyticJac(baVecFunc, x, baAnalyticJ)) return -1;
    time_t baAnalyticEnd = clock();
    std::cout<<"BA Analytic took "
             <<Scalar(baAnalyticEnd - baAnalyticBegin) /
               Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";
    time_t baFFDBegin = clock();
    if (0 != baFFDJac(baVecFunc, x, baFFDJ)) return -1;
    time_t baFFDEnd = clock();
    std::cout<<"BA FFD took "
             <<Scalar(baFFDEnd - baFFDBegin) /
               Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";

    if (baAnalyticJ.m_camNum != baFFDJ.m_camNum ||
        baAnalyticJ.m_ptNum != baFFDJ.m_ptNum ||
        baAnalyticJ.m_camsDrv.size() != baFFDJ.m_camsDrv.size() ||
        baAnalyticJ.m_ptsDrv.size() != baFFDJ.m_ptsDrv.size())
    {
      std::cout<<"Number of cameras or number of points not match.\n";
      return -1;
    }

    Err rst = 0;
    const Scalar threshold = Scalar(1e-5);
    for (Index i = 0; i < Index(baFFDJ.m_camsDrv.size()); i++)
    {
      const typename BAAnaliticJac::CamDrvBlk& baAnalyticCamDrvBlk =
        baAnalyticJ.m_camsDrv[i];
      const typename BAFFDJac::CamDrvBlk& baFFDCamDrvBlk =
        baFFDJ.m_camsDrv[i];

      if (baAnalyticCamDrvBlk.m_drvBlk.rows() != 
          baFFDCamDrvBlk.m_drvBlk.rows() ||
          baAnalyticCamDrvBlk.m_drvBlk.cols() != 
          baFFDCamDrvBlk.m_drvBlk.cols())
      {
        std::cout<<"Size of camera derivative matrix not match.\n";
        rst = -1;
        break;
      }

      for (Index j = 0; j < baFFDCamDrvBlk.m_drvBlk.rows(); j++)
      {
        for (Index k = 0; k < baFFDCamDrvBlk.m_drvBlk.cols(); k++)
        {
          Scalar baFFDValue = baFFDCamDrvBlk.m_drvBlk(j, k);
          Scalar baAnalyticValue = baAnalyticCamDrvBlk.m_drvBlk(j, k);
          Scalar absErr = std::abs(baFFDValue - baAnalyticValue);
          if (absErr > threshold)
          {
            std::cout<<"Deference between the camera derivative matrix value is too big\n";
            std::cout<<"Camera derivative matrix Id:"<<i<<" .\n";
            std::cout<<"baFFDCamDrvBlk["<<j<<", "<<k<<"] = "
                     <<baFFDValue<<".\n";
            std::cout<<"baAnalyticCamDrvBlk["<<j<<", "<<k<<"] = "
                     <<baAnalyticValue<<".\n";
            rst = -1;
          }
        }
      }
    }

    for (Index i = 0; i < Index(baFFDJ.m_ptsDrv.size()); i++)
    {
      const typename BAAnaliticJac::PtDrvBlk& baAnalyticPtDrvBlk =
        baAnalyticJ.m_ptsDrv[i];
      const typename BAFFDJac::PtDrvBlk& baFFDPtDrvBlk =
        baFFDJ.m_ptsDrv[i];

      if (baAnalyticPtDrvBlk.m_drvBlk.rows() != 
          baFFDPtDrvBlk.m_drvBlk.rows() ||
          baAnalyticPtDrvBlk.m_drvBlk.cols() != 
          baFFDPtDrvBlk.m_drvBlk.cols())
      {
        std::cout<<"Size of point derivative matrix not match.\n";
        rst = -1;
        break;
      }

      for (Index j = 0; j < baFFDPtDrvBlk.m_drvBlk.rows(); j++)
      {
        for (Index k = 0; k < baFFDPtDrvBlk.m_drvBlk.cols(); k++)
        {
          Scalar baFFDValue = baFFDPtDrvBlk.m_drvBlk(j, k);
          Scalar baAnalyticValue = baAnalyticPtDrvBlk.m_drvBlk(j, k);
          Scalar absErr = std::abs(baFFDValue - baAnalyticValue);
          if (absErr > threshold)
          {
            std::cout<<"Deference between the point derivative matrix value is too big\n";
            std::cout<<"Point derivative matrix Id:"<<i<<" .\n";
            std::cout<<"baFFDPtDrvBlk["<<j<<", "<<k<<"] = "
                     <<baFFDValue<<".\n";
            std::cout<<"baAnalyticPtDrvBlk["<<j<<", "<<k<<"] = "
                     <<baAnalyticValue<<".\n";
            rst = -1;
          }
        }
      }
    }

    if (withFFD)
    {
      for (Index i = 0; i < ySize; i++)
      {
        for (Index j = 0; j < xSize; j++)
        {
          Scalar ffdS = ffdJ.coeff(i, j);
          Scalar baFFDS = baFFDJ.coeff(i, j);
          Scalar absErr = abs(baFFDS - ffdS);
          if (absErr > Scalar(1e-8))
          {
            std::cout<<"baFFDJ[i, j]="<<baFFDS<<".\n";
            std::cout<<"ffdJ[i, j]="<<ffdS<<".\n";
            return -1;
          }
        }
      }
    }

    return rst;
  }
};

TEST(TestBANaiveJac, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImgDim;

  typedef hs::sfm::ba::BANaiveSyntheticDataGenerator<Scalar, ImgDim> DataGen;
  typedef TestBANaiveJac<Scalar> Test;
  typedef Test::BAVecFunc BAVecFunc;
  typedef BAVecFunc::XVec XVec;
  typedef BAVecFunc::YVec YVec;

  Scalar f = 0.006;
  size_t stripNum = 5;
  size_t camsNumInStrip = 6;
  Scalar grdRes = 0.1;
  ImgDim imgW = 3648;
  ImgDim imgH = 2736;
  Scalar pixSize = 0.00000203311408298266;
  size_t ptsNum = 100;
  Scalar lateralOverlap = 0.6;
  Scalar longitudinalOverlap = 0.8;
  Scalar sceneMaxHeight = 50;
  Scalar camHeightDev = 2;
  Scalar camPlannarDev = 2;
  Scalar camRotDev = 10;
  Scalar nwAngle = 60;

  BAVecFunc baVecFunc;
  XVec x;
  YVec y;

  DataGen dataGen(f, stripNum, camsNumInStrip, grdRes,
          imgW, imgH, pixSize, ptsNum,
          lateralOverlap, longitudinalOverlap,
          sceneMaxHeight,
          camHeightDev, camPlannarDev, camRotDev, nwAngle);
  ASSERT_EQ(0, dataGen(baVecFunc, x, y));

  //ASSERT_EQ(0, Test::test(baVecFunc, x, true));

  f = 0.006;
  stripNum = 1;
  camsNumInStrip = 4;
  grdRes = 0.1;
  imgW = 6000;
  imgH = 4000;
  pixSize = 0.00000203311408298266;
  ptsNum = 10;
  lateralOverlap = 0.99;
  longitudinalOverlap = 0.99;
  sceneMaxHeight = 50;
  camHeightDev = 0.01;
  camPlannarDev = 0.01;
  camRotDev = 0.01;
  nwAngle = 60;

  DataGen dataGen1(f, stripNum, camsNumInStrip, grdRes,
    imgW, imgH, pixSize, ptsNum,
    lateralOverlap, longitudinalOverlap,
    sceneMaxHeight,
    camHeightDev, camPlannarDev, camRotDev, nwAngle);
  ASSERT_EQ(0, dataGen1(baVecFunc, x, y));

  ASSERT_EQ(0, Test::test(baVecFunc, x, true));
}

TEST(TestBANaiveJac, BigDataTest)
{
  typedef double Scalar;
  typedef size_t ImgDim;

  typedef hs::sfm::ba::BANaiveSyntheticDataGenerator<Scalar, ImgDim> DataGen;
  typedef TestBANaiveJac<Scalar> Test;
  typedef Test::BAVecFunc BAVecFunc;
  typedef BAVecFunc::XVec XVec;
  typedef BAVecFunc::YVec YVec;

  Scalar f = 0.006;
  size_t stripNum = 5;
  size_t camsNumInStrip = 40;
  Scalar grdRes = 0.1;
  ImgDim imgW = 3648;
  ImgDim imgH = 2736;
  Scalar pixSize = 0.00000203311408298266;
  size_t ptsNum = 100000;
  Scalar lateralOverlap = 0.6;
  Scalar longitudinalOverlap = 0.8;
  Scalar sceneMaxHeight = 50;
  Scalar camHeightDev = 2;
  Scalar camPlannarDev = 2;
  Scalar camRotDev = 10;
  Scalar nwAngle = 60;

  BAVecFunc baVecFunc;
  XVec x;
  YVec y;

  DataGen dataGen(f, stripNum, camsNumInStrip, grdRes,
          imgW, imgH, pixSize, ptsNum,
          lateralOverlap, longitudinalOverlap,
          sceneMaxHeight,
          camHeightDev, camPlannarDev, camRotDev, nwAngle);
  ASSERT_EQ(0, dataGen(baVecFunc, x, y));

  ASSERT_EQ(0, Test::test(baVecFunc, x, false));

}

}
