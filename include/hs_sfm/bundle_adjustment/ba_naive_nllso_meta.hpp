#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_NLLSO_META_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_NLLSO_META_HPP_

#include "hs_math/linear_algebra/lafunc/arithmetic.hpp"
#include "hs_math/linear_algebra/lafunc/arithmetic_eigen.hpp"
#include "hs_math/linear_algebra/latraits/mat_eigen.hpp"
#include "hs_math/linear_algebra/latraits/vec_eigen.hpp"
#include "hs_optimizor/nllso/meta_fwd.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_analytic_jac.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_feat_cov_inv.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_jac_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_nllso_augmentor.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_nllso_normal_equation.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_nllso_normal_equation_builder.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_nllso_normal_equation_solver.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_residuals.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_vec_func.hpp"

namespace hs
{
namespace optimizor
{
namespace nllso
{

template <typename _Scalar>
struct JacobiType<hs::sfm::ba::BANaiveVecFunc<_Scalar> >
{
  typedef hs::sfm::ba::BANaiveAnalyticJacobian<
            hs::sfm::ba::BANaiveVecFunc<_Scalar> > type;
};

template <typename _Scalar>
struct ResidualType<hs::sfm::ba::BANaiveVecFunc<_Scalar> >
{
  typedef hs::sfm::ba::BANaiveResidualsCalc<
            hs::sfm::ba::BANaiveVecFunc<_Scalar> > type;
};

template <typename _Scalar, typename _Index,
          _Index params_per_cam,
          _Index params_per_pt>
struct NormalEqSolverType<hs::sfm::ba::BANaiveAugmentNormalMatrix<
                            _Scalar, _Index,
                            params_per_cam,
                            params_per_pt>,
                          hs::sfm::ba::BANaiveBVec<
                            _Scalar, _Index,
                            params_per_cam,
                            params_per_pt> >
{
  typedef hs::sfm::ba::BANaiveNormalEquationSolver<_Scalar> type;
};

template <typename _Scalar, typename _Index,
          _Index params_per_feat,
          _Index params_per_cam,
          _Index params_per_pt>
struct NormalEqBuilderType<EIGEN_VEC(_Scalar, Eigen::Dynamic),
                           EIGEN_VEC(_Scalar, Eigen::Dynamic),
                           hs::sfm::ba::BANaiveJacMatrix<
                             _Scalar, _Index,
                             params_per_feat,
                             params_per_cam,
                             params_per_pt>,
                           hs::sfm::ba::BANaiveFeatCovInv<_Scalar> >
{
  typedef hs::sfm::ba::BANaiveNormalEquationBuilder<_Scalar> type;
};

template <typename _Scalar, typename _Index,
          _Index params_per_cam,
          _Index params_per_pt>
struct NormalAugmentorType<hs::sfm::ba::BANaiveNormalMatrix<
                             _Scalar, _Index,
                             params_per_cam,
                             params_per_pt> >
{
  typedef hs::sfm::ba::BANaiveAugmentor<_Scalar> type;
};

}
}

}

#endif
