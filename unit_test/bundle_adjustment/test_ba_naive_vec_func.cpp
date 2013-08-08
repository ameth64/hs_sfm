#include <gtest/gtest.h>

#include "test_ba_naive_base.hpp"

namespace
{

TEST(TestBANaiveVecFunc, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImgDim;

  typedef SyntheticDataGenerator<Scalar, ImgDim> DataGen;

  typedef hs::sfm::ba::BANaiveVecFunc<Scalar> BAVecFunc;
  typedef BAVecFunc::Index Index;
  typedef BAVecFunc::XVec XVec;
  typedef BAVecFunc::YVec YVec;

  Scalar f = 0.006;
  size_t stripNum = 4;
  size_t camsNumInStrip = 5;
  Scalar grdRes = 0.1;
  ImgDim imgW = 3648;
  ImgDim imgH = 2736;
  Scalar pixSize = 0.00000203311408298266;
  size_t ptsNum = 2000;
  Scalar lateralOverlap = 0.6;
  Scalar longitudinalOverlap = 0.8;
  Scalar sceneMaxHeight = 50;
  Scalar camHeightDev = 2;
  Scalar camPlannarDev = 2;
  Scalar camRotDev = 10;
  Scalar nwAngle = 60;

  BAVecFunc baVecFunc;
  XVec x;
  YVec y, y1;

  DataGen dataGen(f, stripNum, camsNumInStrip, grdRes, 
          imgW, imgH, pixSize, ptsNum,
          lateralOverlap, longitudinalOverlap,
          sceneMaxHeight,
          camHeightDev, camPlannarDev, camRotDev, nwAngle);
  ASSERT_EQ(0, dataGen(baVecFunc, x, y));

  ASSERT_EQ(0, baVecFunc(x, y1));

  for (Index i = 0; i < y1.rows(); i++)
  {
    ASSERT_NEAR(y[i], y1[i], Scalar(1e-5));
  }
}

}
