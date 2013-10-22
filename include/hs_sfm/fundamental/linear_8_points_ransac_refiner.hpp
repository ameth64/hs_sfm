#ifndef _HS_SFM_FUNDAMENTAL_LINEAR_8_POINTS_RANSAC_REFINER_HPP_
#define _HS_SFM_FUNDAMENTAL_LINEAR_8_POINTS_RANSAC_REFINER_HPP_

#include <utility>

#include "hs_fit/ransac/ransac.hpp"

#include "hs_sfm/fundamental/linear_8_points_calculator.hpp"

namespace hs
{

namespace fit
{

template <typename _Scalar>
struct PointSetType<std::pair<EIGEN_VECTOR(_Scalar, 2),
                              EIGEN_VECTOR(_Scalar, 2)> >
{
  typedef std::pair<EIGEN_VECTOR(_Scalar, 2),
                    EIGEN_VECTOR(_Scalar, 2)> Point;
  typedef EIGEN_STD_VECTOR(Point) type;
};

}

namespace sfm
{
namespace fundamental
{

template <typename _Scalar>
class Linear8PointsRansacRefiner
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef std::pair<Key, Key> KeyPair;
  typedef EIGEN_STD_VECTOR(KeyPair) KeyPairContainer;

private:
  typedef Linear8PointsCalculator<Scalar> Calculator;
  typedef typename Calculator::FMatrix FMatrix;
  typedef typename Calculator::KeyContainer KeyContainer;

  class FMatrixCalculator
  {
  public:
    typedef FMatrix RansacModel;
    static const size_t model_min_set_size_ = 8;

  public:
    int operator() (const KeyPairContainer& key_pairs, FMatrix& f_matrix) const
    {
      size_t number_of_keys = key_pairs.size();
      if (number_of_keys < model_min_set_size_)
      {
        return -1;
      }

      KeyContainer keys_left(number_of_keys);
      KeyContainer keys_right(number_of_keys);
      for (size_t i = 0; i < number_of_keys; i++)
      {
        keys_left[i] = key_pairs[i].first;
        keys_right[i] = key_pairs[i].second;
      }

      Calculator calculator;
      return calculator(keys_left, keys_right, f_matrix);
    }
  };

  class FMatrixResidualsCalculator
  {
  public:
    typedef Scalar Distance;

  private:
    typedef EIGEN_VECTOR(Scalar, 3) HKey;
    typedef EIGEN_VECTOR(Scalar, 3) HLine;

  public:
    int operator() (const KeyPair& key_pair, const FMatrix& f_matrix,
                    Distance& distance) const
    {
      HKey hkey_left;
      hkey_left << key_pair.first[0],
                   key_pair.first[1],
                   Scalar(1);
      HKey hkey_right;
      hkey_right << key_pair.second[0],
                    key_pair.second[1],
                    Scalar(1);

      HLine epiline_left = f_matrix * hkey_left;
      HLine epiline_right = f_matrix.transpose() * hkey_right;
      Scalar error_left = std::abs(hkey_right.dot(epiline_left)) /
                          epiline_left.segment(0, 2).norm();
      Scalar error_right = std::abs(hkey_left.dot(epiline_right)) /
                           epiline_right.segment(0, 2).norm();

      distance = (error_left + error_right) / Scalar(2);

      return 0;
    }
  };

  typedef hs::fit::Ransac<KeyPair,
                          FMatrixCalculator,
                          FMatrixResidualsCalculator> RansacRefiner;
  
public:
  typedef typename RansacRefiner::IndexSet IndexSet;

public:
  Err operator() (const KeyPairContainer& key_pairs,
                  Scalar distance_threshold,
                  KeyPairContainer& refined_key_pairs,
                  IndexSet& inlier_indices) const
  {
    FMatrixCalculator model_calculator;
    FMatrixResidualsCalculator distance_calculator;
    RansacRefiner ransac_refiner(model_calculator, distance_calculator);
    
    ransac_refiner.SetAlphaThreshold(Scalar(0.95));
    ransac_refiner.SetDistanceThreshold(distance_threshold);
    KeyPairContainer best_key_pairs;
    return (ransac_refiner(key_pairs, refined_key_pairs,
                           inlier_indices, best_key_pairs));
  }
};

}
}
}

#endif