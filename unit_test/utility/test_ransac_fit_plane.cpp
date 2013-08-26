#include <gtest/gtest.h>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_sfm/utility/ransac_fit_plane.hpp"

namespace
{

template <typename _Scalar>
class TestRansacFitPlannar
{
public:
  typedef _Scalar Scalar;
  typedef EIGEN_VEC(Scalar, 3) Pt;
  EIGEN_BASE_DEF(Scalar);

  /**
   *  测试ransac拟合三维点平面。
   *  平面真值为xy平面作相似变换后的平面。
   */
  static int testRansacFitPlannar(const Mat33& rot,
               const Vec3& trans,
               Scalar scale,
               const Vec2& min,
               const Vec2& max,
               Scalar sigma,
               Scalar outlierRatio,
               Scalar outlierThd,
               size_t ptsNum)
  {
    Mat33 cov = Mat33::Identity();
    cov *= sigma;
    std::vector<Pt> pts;
    for (size_t i = 0; i < ptsNum; i++)
    {
      Vec2 plannarPt;
      hs::math::random::UniformRandomVar<Scalar, 2>::uniformRandomVar(min, max, 
                              plannarPt);

      Vec3 mean;
      mean << plannarPt[0],
          plannarPt[1],
          0;
      Vec3 pt;
      hs::math::random::NormalRandomVar<Scalar, 3>::normRandomVar(mean, cov, pt);

      Scalar rnd;
      hs::math::random::UniformRandomVar<Scalar, 1>::uniformRandomVar(Scalar(0),
                                Scalar(1),
                                rnd);
      if (rnd < outlierRatio)
      {
        Scalar outlier;
        hs::math::random::UniformRandomVar<Scalar, 1>::uniformRandomVar(
          outlierThd, outlierThd * 3, outlier);

        pt[0] += outlier;
        pt[1] += outlier;
        pt[2] += outlier;
      }

      pt = rot * pt * scale + trans;

      Pt p;
      p[0] = pt[0];
      p[1] = pt[1];
      p[2] = pt[2];
      pts.push_back(p);
    }

    Vec4 pln;
    hs::sfm::RansacFitPlane<Pt> ransacFit;
    if (ransacFit(pts, pln, outlierThd) != 0) return -1;

    size_t outlierNum = 0;
    for (size_t i = 0; i < ptsNum; i++)
    {
      Vec4 hPt;
      hPt << pts[i][0],
           pts[i][1],
           pts[i][2],
           1;
      Scalar dist = hPt.dot(pln) / pln.segment(0, 3).norm();
      if (dist > outlierThd * scale)
      {
        outlierNum++;
      }
    }

    return Scalar(outlierNum) / Scalar(ptsNum) < outlierRatio * 2 ? 0 : -1;
  }
};

TEST(TestRansacFitPlannar, SimpleTest)
{
  typedef double Scalar;
  EIGEN_BASE_DEF(Scalar);

  Mat33 rot = Mat33::Identity();
  Vec3 trans;
  trans << 0,
       0,
       10;
  Scalar scale = 5.0;

  Vec2 min;
  min << -10,
       -10;
  Vec2 max;
  max << 10,
       10;

  Scalar sigma = 1.0;
  Scalar outlierRatio = 0.2;
  Scalar outlierThd = 5;
  size_t ptsNum = 100;

  ASSERT_EQ(0, TestRansacFitPlannar<Scalar>::
    testRansacFitPlannar(rot, 
               trans,
               scale,
               min,
               max,
               sigma,
               outlierRatio,
               outlierThd,
               ptsNum));
}

}
