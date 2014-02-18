#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_MATRIX_MAX_DIAGONAL_VALUE_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_MATRIX_MAX_DIAGONAL_VALUE_CALCULATOR_HPP_

#include <algorithm>
#include <limits>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BANaiveNormalMatrixMaxDiagonalValueCalculator
{
public:
  typedef _Scalar Scalar;
  typedef BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef BANaiveNormalMatrix<Scalar, Index,
                              VectorFunction::params_per_camera_,
                              VectorFunction::params_per_point_>
          NormalMatrix;

  Scalar operator() (const NormalMatrix& normal_matrix) const
  {
    Index rows = normal_matrix.number_of_cameras *
                 VectorFunction::params_per_camera_ +
                 normal_matrix.number_of_points *
                 VectorFunction::params_per_point_;
    Scalar max_value = -std::numeric_limits<Scalar>::max();
    for (Index i = 0; i < rows; i++)
    {
      max_value = std::max(max_value, normal_matrix.coeff(i, i));
    }

    return max_value;
  }
};

}
}
}

#endif
