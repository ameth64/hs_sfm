#ifndef _HS_SFM_ESSENTIAL_EMATRIX_5_POINTS_RANSAC_REFINER_HPP_
#define _HS_SFM_ESSENTIAL_EMATRIX_5_POINTS_RANSAC_REFINER_HPP_

#include <utility>

#include "hs_fit/ransac/ransac.hpp"

#include "hs_sfm/essential/ematrix_5_points_calculator.hpp"

namespace hs
{
namespace fit
{

template <typename _Scalar>
struct PointSetType<std::pair<EIGEN_VECTOR(_Scalar, 3),
                              EIGEN_VECTOR(_Scalar, 3)> >
{
  typedef std::pair<EIGEN_VECTOR(_Scalar, 3),
                    EIGEN_VECTOR(_Scalar, 3)> Point;
  typedef EIGEN_STD_VECTOR(Point) type;
};

}

namespace sfm
{
namespace essential
{

template <typename _Scalar>
class EMatrix5PointsRansacRefiner
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  
private:
  typedef EMatrix5PointsCalculator<Scalar> Calculator;
public:
  typedef typename Calculator::HKey HKey;
  typedef typename Calculator::HKeyPair HKeyPair;
  typedef typename Calculator::HKeyPairContainer HKeyPairContainer;
private:
  typedef typename Calculator::EMatrixHypotheses EMatrixHypotheses;

  class RansacModelCalculator
  {
  public:
    typedef EMatrixHypotheses RansacModel;
    static const size_t model_min_set_size_ = 5;

  private:
    typedef typename Calculator::HypothesesGenerator Generator;

  public:
    int operator() (const HKeyPairContainer& key_pairs,
                    EMatrixHypotheses& ematrix_hypotheses) const
    {
      Generator generator;
      return generator(key_pairs, ematrix_hypotheses);

      return 0;
    }
  };

  class RansacDistanceCalculator
  {
  public:
    typedef Scalar Distance;

  private:
    typedef typename Calculator::HypothesesEvaluator Evaluator;

  public:
    int operator() (const HKeyPair& key_pair,
                    const EMatrixHypotheses& ematrix_hypotheses,
                    Distance& distance) const
    {
      Evaluator evaluator;
      size_t best_hypothesis;
      return evaluator(key_pair, ematrix_hypotheses,
                       best_hypothesis, distance);
    }
  };

  typedef hs::fit::Ransac<HKeyPair,
                          RansacModelCalculator,
                          RansacDistanceCalculator> RansacRefiner;

public:
  typedef typename RansacRefiner::IndexSet IndexSet;

public:
  Err operator() (const HKeyPairContainer& key_pairs,
                  Scalar distance_threshold,
                  HKeyPairContainer& refined_key_pairs,
                  IndexSet& inlier_indices) const
  {
    RansacModelCalculator model_calculator;
    RansacDistanceCalculator distance_calculator;
    RansacRefiner ransac_refiner(model_calculator, distance_calculator);

    ransac_refiner.SetAlphaThreshold(Scalar(0.95));
    ransac_refiner.SetDistanceThreshold(distance_threshold);
    HKeyPairContainer best_key_pairs;
    return (ransac_refiner(key_pairs, refined_key_pairs,
                           inlier_indices, best_key_pairs));
  }

};

}
}
}

#endif