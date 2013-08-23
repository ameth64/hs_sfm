#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_XCOV_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_XCOV_CALCULATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_nllso_normal_equation.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VecFunc>
class BANaiveXCovCalculator;

template <typename _Scalar>
class BANaiveXCovCalculator<BANaiveVecFunc<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef BANaiveVecFunc<Scalar> VecFunc;
  typedef typename VecFunc::Index Index;
  typedef typename VecFunc::XVec XVec;
  typedef typename VecFunc::YVec YVec;

  typedef BANaiveNormalMatrix<Scalar, Index,
                              VecFunc::m_paramsPerCam,
                              VecFunc::m_paramsPerPt> NormalMatrix;

  enum
  {
    EigenDefaultMajor =
#if defined(EIGEN_DEFAULT_TO_ROW_MAJOR)
    Eigen::RowMajor
#else
    Eigen::ColMajor
#endif
  };
  typedef EIGEN_MAT(Scalar, Eigen::Dynamic, Eigen::Dynamic) DenseNormalMatrix;
  typedef DenseNormalMatrix XCovariance;

  typedef int Err;

  Err operator()(const NormalMatrix& normal_matrix,
                 XCovariance& x_covariance) const
  {
    Index camera_params_size = normal_matrix.number_of_cameras *
                               VecFunc::m_paramsPerCam;
    Index point_params_size = normal_matrix.number_of_points *
                              VecFunc::m_paramsPerPt;
    Index x_size = camera_params_size + point_params_size;

    DenseNormalMatrix dense_normal_matrix(x_size, x_size);
    for (Index i = 0; i < x_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        dense_normal_matrix(i, j) = normal_matrix.coeff(i, j);
      }
    }

    //计算正规矩阵的伪逆即为X的协方差矩阵
    Eigen::JacobiSVD<DenseNormalMatrix> svd(dense_normal_matrix);
    DenseNormalMatrix U = svd.matrixU();
    DenseNormalMatrix V = svd.matrixV();
    DenseNormalMatrix D = DenseNormalMatrix::Zero(x_size, x_size);
    EIGEN_VEC(Scalar, Eigen::Dynamic) d = svd.singularValues();
    const Scalar precision = Scalar(1e-6);
    for (Index i = 0; i < x_size; i++)
    {
      if (std::abs(d(i)) > precision)
      {
        D(i, i) = Scalar(1) / d(i);
      }
    }

    x_covariance.resize(x_size, x_size);
    x_covariance = V * D * U.transpose();

    return 0;
  }
};

}
}
}

#endif