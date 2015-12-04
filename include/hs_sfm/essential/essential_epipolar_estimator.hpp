#ifndef _HS_SFM_ESSENTIAL_ESSENTIAL_EPIPOLAR_ESTIMATOR_HPP_
#define _HS_SFM_ESSENTIAL_ESSENTIAL_EPIPOLAR_ESTIMATOR_HPP_

#include "hs_sfm/essential/ematrix_5_points_calculator.hpp"
#include "hs_sfm/essential/ematrix_5_points_ransac_refiner.hpp"
#include "hs_sfm/essential/ematrix_extrinsic_params_points_calculator.hpp"

namespace hs
{
namespace sfm
{
namespace essential
{

template <typename _Scalar>
class EssentialEpipolarEstimator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef EMatrix5PointsCalculator<Scalar> EMatrixCalculator;
  typedef EMatrix5PointsRansacRefiner<Scalar> EMatrixRefiner;
  typedef EMatrixExtrinsicParamsPointsCalculator<Scalar>
          EMatrixStructureCalculator;
  typedef typename EMatrixCalculator::HKey HKey;
  typedef typename EMatrixCalculator::HKeyPair HKeyPair;
  typedef typename EMatrixCalculator::HKeyPairContainer HKeyPairContainer;
  typedef typename EMatrixStructureCalculator::ExtrinsicParams
          ExtrinsicParams;
  typedef typename EMatrixStructureCalculator::Point Point;
  typedef typename EMatrixStructureCalculator::PointContainer PointContainer;

public:
  Err EstimateRobust(const HKeyPairContainer& hkey_pairs,
                     ExtrinsicParams& extrinsic_params,
                     PointContainer& points) const
  {
    return 0;
  }
};

}
}
}

#endif
