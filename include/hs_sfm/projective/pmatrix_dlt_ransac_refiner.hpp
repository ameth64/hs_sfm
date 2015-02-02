#ifndef _HS_SFM_PROJECTIVE_PMATRIX_DLT_RANSAC_REFINER_HPP_
#define _HS_SFM_PROJECTIVE_PMATRIX_DLT_RANSAC_REFINER_HPP_

#include <utility>

#include "hs_fit/ransac/ransac.hpp"

#include "hs_sfm/projective/pmatrix_dlt_calculator.hpp"

namespace hs
{

namespace fit
{

template <typename _Scalar>
struct PointSetType<std::pair<EIGEN_VECTOR(_Scalar, 2),
                              EIGEN_VECTOR(_Scalar, 3)> >
{
  typedef std::pair<EIGEN_VECTOR(_Scalar, 2),
                    EIGEN_VECTOR(_Scalar, 3)> Point;
  typedef EIGEN_STD_VECTOR(Point) type;
};

}

namespace sfm
{
namespace projective
{

template <typename _Scalar>
class PMatrixDLTRansacRefiner
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  
private:
  typedef PMatrixDLTCalculator<Scalar> Calculator;
public:
  typedef typename Calculator::PMatrix PMatrix;
  typedef typename Calculator::Correspondence Correspondence;
  typedef typename Calculator::CorrespondenceContainer CorrespondenceContainer;

private:
  class RansacModelCalculator
  {
  public:
    typedef PMatrix RansacModel;
    static const size_t model_min_set_size_ = 8;

  public:
    int operator() (const CorrespondenceContainer& correspondences,
                    PMatrix& p_matrix) const
    {
      Calculator calculator;
      return calculator(correspondences, p_matrix);
    }
  };

  class RansacDistanceCalculator
  {
  public:
    typedef Scalar Distance;

  private:
    typedef EIGEN_VECTOR(Scalar, 3) HKey;
    typedef EIGEN_VECTOR(Scalar, 4) HPoint;

  public:
    int operator()(const Correspondence& correspondence,
                   const PMatrix& p_matrix,
                   Distance& distance) const
    {
      HPoint hpoint;
      hpoint.template segment<3>(0) = correspondence.second;
      hpoint[3] = Scalar(1);

      HKey hkey = p_matrix * hpoint;
      hkey /= hkey(2);

      distance = (hkey.template segment<2>(0) -
                  correspondence.first).norm();

      //std::cout<<"distance:"<<distance<<"\n";
      return 0;
    }
  };

  typedef hs::fit::Ransac<Correspondence,
                          RansacModelCalculator,
                          RansacDistanceCalculator> RansacRefiner;

public:
  typedef typename RansacRefiner::IndexSet IndexSet;

public:
  Err operator() (const CorrespondenceContainer& correspondences,
                  Scalar distance_threshold,
                  CorrespondenceContainer& refined_correspondences,
                  IndexSet& inlier_indices) const
  {
    RansacModelCalculator model_calculator;
    RansacDistanceCalculator distance_calculator;
    RansacRefiner refiner(model_calculator, distance_calculator);

    refiner.SetAlphaThreshold(Scalar(0.95));
    refiner.SetDistanceThreshold(distance_threshold);
    CorrespondenceContainer best_correspondences;
    return (refiner(correspondences, refined_correspondences,
                    inlier_indices, best_correspondences));
  }
};

}
}

}

#endif
