#include <iostream>
#include <ctime>

#include <gtest/gtest.h>

#include "hs_math/linear_algebra/latraits/vec_eigen.hpp"
#include "hs_math/linear_algebra/lafunc/arithmetic_eigen.hpp"
#include "hs_math/fdjac/ffd_sparse_jac.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_jac.hpp"

#include "test_ba_naive_base.hpp"

namespace
{
  TEST(TestBANaiveJac, SimpleTest)
  {
    typedef double Scalar;
    typedef size_t ImgDim;

    typedef SyntheticDataGenerator<Scalar, ImgDim> DataGen;

    typedef hs::sfm::ba::BANaiveVecFunc<Scalar> BAVecFunc;
    typedef BAVecFunc::Index Index;
    typedef BAVecFunc::XVec XVec;
    typedef BAVecFunc::YVec YVec;

    typedef hs::math::fdjac::FwdFiniteDiffSparseJacobian<BAVecFunc> FFDJacobian;
    typedef FFDJacobian::Jac FFDJac;
    typedef hs::sfm::ba::BANaiveJacobian<BAVecFunc> BAJacobian;
    typedef BAJacobian::Jac BAJac;

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

    FFDJacobian ffdJac;
    BAJacobian baJac;

    FFDJac ffdJ;
    BAJac baJ;

    time_t ffdBegin = clock();
    ASSERT_EQ(0, ffdJac(baVecFunc, x, ffdJ));
    time_t ffdEnd = clock();
    time_t baBegin = clock();
    ASSERT_EQ(0, baJac(baVecFunc, x, baJ));
    time_t baEnd = clock();

    std::cout<<"FFD took "<<Scalar(ffdEnd - ffdBegin) / 
                Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";
    std::cout<<"BA took "<<Scalar(baEnd - baBegin) /
                 Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";

    for (Index i = 0; i < ffdJ.rows(); i++)
    {
      for (Index j = 0; j < ffdJ.cols(); j++)
      {
        Scalar ffdS = ffdJ.coeff(i, j);
        Scalar baS = baJ.coeff(i, j);

        Scalar absErr = abs(baS - ffdS);
        ASSERT_NEAR(absErr, 0, 1e-4);
      }
    }

    int bp = 0;
  }
}
