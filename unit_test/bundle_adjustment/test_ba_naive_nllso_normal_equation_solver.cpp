#include <iostream>

#include <gtest/gtest.h>

#include "hs_sfm/bundle_adjustment/ba_naive_analytic_jac.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_nllso_normal_equation_builder.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_synthetic_data_generator.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_nllso_augmentor.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_nllso_normal_equation_solver.hpp"

namespace
{

template <typename _Scalar>
class TestBANaiveNLLSONormalEquationSolver
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef hs::sfm::ba::BANaiveVecFunc<Scalar> BAVecFunc;
  typedef typename BAVecFunc::XVec XVec;
  typedef typename BAVecFunc::YVec YVec;
  typedef typename BAVecFunc::Index Index;
  typedef hs::sfm::ba::BANaiveNormalEquationBuilder<Scalar>
          NormalEquationBuilder;
  typedef typename NormalEquationBuilder::YCovInv YCovInv;
  typedef typename NormalEquationBuilder::NormalMat NormalMat;
  typedef typename NormalEquationBuilder::BVec BVec;
  typedef typename NormalEquationBuilder::ResidualsCalc ResidualsCalc;
  typedef typename ResidualsCalc::Residuals Residuals;

  typedef hs::sfm::ba::BANaiveAnalyticJacobian<BAVecFunc> BAAnalyticJacobian;
  typedef typename BAAnalyticJacobian::Jac BAAnalyticJac;

  typedef hs::sfm::ba::BANaiveNormalEquationSolver<Scalar>
          NormalEquationSolver;
  typedef typename NormalEquationSolver::AugmentNormalMat AugmentNormalMat;

  typedef hs::sfm::ba::BANaiveAugmentor<Scalar> Augmentor;
  
  typedef EIGEN_MAT(Scalar, Eigen::Dynamic, Eigen::Dynamic)
          MatXX;
  typedef EIGEN_MAT(Scalar, 2, 2) Mat22;

  static Err Test(const BAVecFunc& f, const XVec& x, const Mat22& feat_cov)
  {
    Err result = 0;

    YVec y;
    if (f(x, y) != 0)
    {
      std::cout<<"Fail to call vector function.\n";
      return -1;
    }

    Index feat_num = f.getFeatNum();
    MatXX dense_y_cov, dense_y_cov_inv;
    YCovInv y_cov_inv;
    if (GenFeatCovInv(feat_cov, feat_num,
      dense_y_cov, dense_y_cov_inv,
      y_cov_inv) != 0)
    {
      std::cout<<"GenFeatCovInv Failed.\n";
      return -1;
    }

    YVec noised_y = y;
    if (GenNoisedYVec(dense_y_cov, noised_y) != 0)
    {
      std::cout<<"GenNoisedYVec Failed.\n";
      return -1;
    }

    BAAnalyticJacobian ba_analytic_jac;
    BAAnalyticJac ba_analytic_j;
    if (ba_analytic_jac(f, x, ba_analytic_j) != 0)
    {
      std::cout<<"ba_analytic_jac failed.\n";
      return -1;
    }

    ResidualsCalc residuals_calc;
    Residuals r = residuals_calc(noised_y, y);

    NormalMat N;
    BVec b;
    NormalEquationBuilder builder;
    if (builder(ba_analytic_j, r, y_cov_inv, N, b) != 0)
    {
      std::cout<<"builder failed.\n";
      return -1;
    }

    Scalar max_diag = hs::math::la::MatMaxDiagValFunc<NormalMat>()(N);

    Scalar mu = /*Scalar(1e-3) * */max_diag;

    Augmentor augmentor;
    AugmentNormalMat AN = augmentor(N, mu);

    NormalEquationSolver solver;
    XVec x_delta;
    if (solver(AN, b, x_delta) != 0)
    {
      std::cout<<"solver failed.\n";
      return -1;
    }

    MatXX dense_N;
    XVec dense_b;
    Index x_size = f.getXSize();
    dense_N.resize(x_size, x_size);
    dense_b.resize(x_size);
    for (Index i = 0; i < x_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        dense_N(i, j) = N.coeff(i, j);
      }
      dense_b[i] = b[i];
      dense_N(i, i) += mu;
    }

    Scalar threshold = Scalar(1e-8);
    
    if ((dense_N * x_delta).isApprox(dense_b, threshold))
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }

private:
  static Err GenFeatCovInv(const Mat22& feat_cov,
                           Index feat_num,
                           MatXX& dense_y_cov,
                           MatXX& dense_y_cov_inv,
                           YCovInv& y_cov_inv)
  {
    y_cov_inv.m_blocks.clear();
    dense_y_cov_inv.resize(feat_num * BAVecFunc::m_paramsPerFeat,
                           feat_num * BAVecFunc::m_paramsPerFeat);
    dense_y_cov_inv.setZero();
    dense_y_cov.resize(feat_num * BAVecFunc::m_paramsPerFeat,
                       feat_num * BAVecFunc::m_paramsPerFeat);
    dense_y_cov.setZero();
    for (Index i = 0; i < feat_num; i++)
    {
      dense_y_cov.template block<BAVecFunc::m_paramsPerFeat,
                                 BAVecFunc::m_paramsPerFeat>(
        i * BAVecFunc::m_paramsPerFeat,
        i * BAVecFunc::m_paramsPerFeat) = feat_cov;
      dense_y_cov_inv.template block<BAVecFunc::m_paramsPerFeat,
                                     BAVecFunc::m_paramsPerFeat>(
        i * BAVecFunc::m_paramsPerFeat,
        i * BAVecFunc::m_paramsPerFeat) = feat_cov.inverse();
      y_cov_inv.m_blocks.push_back(feat_cov.inverse());
    }

    return 0;
  }

  static Err GenNoisedYVec(const MatXX& dense_y_cov,
                           YVec& y)
  {
    YVec mean = y;
    return hs::math::random::NormalRandomVar<Scalar, 
                                             Eigen::Dynamic>::normRandomVar(
      mean, dense_y_cov, y);
  }
};

TEST(TestBANaiveNLLSONormalEquationSolver, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImgDim;

  typedef hs::sfm::ba::BANaiveSyntheticDataGenerator<Scalar, ImgDim> DataGen;
  typedef TestBANaiveNLLSONormalEquationSolver<Scalar> Test;
  typedef Test::Mat22 Mat22;
  typedef Test::BAVecFunc BAVecFunc;
  typedef BAVecFunc::XVec XVec;
  typedef BAVecFunc::YVec YVec;

  Scalar f = 0.006;
  size_t strip_num = 1;
  size_t cams_num_in_strip = 4;
  Scalar grd_res = 0.1;
  ImgDim img_w = 6000;
  ImgDim img_h = 4000;
  Scalar pix_size = 0.00000203311408298266;
  size_t pts_num = 10;
  Scalar lateral_overlap = 0.99;
  Scalar longitudinal_overlap = 0.99;
  Scalar scene_max_height = 50;
  Scalar cam_height_dev = 0.01;
  Scalar cam_plannar_dev = 0.01;
  Scalar cam_rot_dev = 0.01;
  Scalar nw_angle = 60;

  BAVecFunc ba_vec_func;
  XVec x;
  YVec y;

  DataGen data_gen(f, strip_num, cams_num_in_strip, grd_res,
    img_w, img_h, pix_size, pts_num,
    lateral_overlap, longitudinal_overlap,
    scene_max_height,
    cam_height_dev, cam_plannar_dev, cam_rot_dev, nw_angle);
  ASSERT_EQ(0, data_gen(ba_vec_func, x, y));
  Scalar f_in_pix = f / pix_size;

  Mat22 feat_cov = Mat22::Identity();
  feat_cov *= Scalar(1) / (f_in_pix * f_in_pix);

  ASSERT_EQ(0, Test::Test(ba_vec_func, x, feat_cov));
}

}
