#ifndef _HS_SFM_ESSENTIAL_EMATRIX_5_POINTS_CALCULATOR_HPP_
#define _HS_SFM_ESSENTIAL_EMATRIX_5_POINTS_CALCULATOR_HPP_

#include <utility>
#include <limits>
#include <vector>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_math/polynomial/multi_polynomial.hpp"

namespace hs
{
namespace sfm
{
namespace essential
{

template <typename _Scalar>
class EMatrix5PointsCalculator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef EIGEN_MATRIX(Scalar, 3, 3) EMatrix;
  typedef EIGEN_STD_VECTOR(EMatrix) EMatrixHypotheses;
  typedef EIGEN_VECTOR(Scalar, 3) HKey;
  typedef std::pair<HKey, HKey> HKeyPair;
  typedef EIGEN_STD_VECTOR(HKeyPair) HKeyPairContainer;

private:
  typedef EIGEN_VECTOR(Scalar, 3) HLine;

public:

  class HypothesesGenerator	//!
  {
  private:
    typedef EIGEN_MATRIX(Scalar, 9, 4) NullspaceBasis;
    typedef EIGEN_MATRIX(Scalar, 10, 10) GroebnerBasis;
    typedef EIGEN_MATRIX(Scalar, 10, 10) ActionMatrix;

    typedef hs::math::MultiPolynomial<Scalar, 3, 3> Polynomial33;

  public:
    Err operator() (const HKeyPairContainer& key_pairs,
                    EMatrixHypotheses& ematrix_hypotheses) const
    {
      if (key_pairs.size() < 5)
      {
        return -1;
      }

      NullspaceBasis basis;
      if (ComputeNullspaceBasis(key_pairs, basis) != 0)
      {
        return -1;
      }

      std::vector<Polynomial33> constraint(10);
      if (ComputeConstraintMatrix(basis, constraint) != 0)
      {
        return -1;
      }

      GroebnerBasis groebner_basis;
      if (ComputeGroebnerBasis(constraint, groebner_basis) != 0)
      {
        return -1;
      }

      ActionMatrix action_matrix;
      if (ComputeActionMatrix(groebner_basis, action_matrix) != 0)
      {
        return -1;
      }

      if (ComputeEMatrixHypothese(action_matrix, basis,
                                  ematrix_hypotheses) != 0)
      {
        return -1;
      }

      return 0;
    }

  private:
    Err ComputeNullspaceBasis(const HKeyPairContainer& key_pairs,
                              NullspaceBasis& basis) const	//求解零空间的一组基
    {
      typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic)
              KeyPairsMatrix;
      size_t number_of_key_pairs = key_pairs.size();
      KeyPairsMatrix key_pairs_matrix(number_of_key_pairs, 9);

      //排列
      for (size_t i = 0; i < number_of_key_pairs; i++)
      {
        //l for left, r for right.
        const HKey& l = key_pairs[i].first;
        const HKey& r = key_pairs[i].second;

        key_pairs_matrix.row(i) << l[0] * r[0], l[1] * r[0], l[2] * r[0],
                                   l[0] * r[1], l[1] * r[1], l[2] * r[1],
                                   l[0] * r[2], l[1] * r[2], l[2] * r[2];
      }

      //计算svd，求出最小的四个奇异值对应的奇异向量
      size_t number_of_v_columns = std::min<size_t>(number_of_key_pairs, 9);
      Eigen::JacobiSVD<KeyPairsMatrix> svd(key_pairs_matrix,
                                           Eigen::ComputeFullV);
      for (size_t i = 0; i < 4; i++)
      {
        basis.col(i) = svd.matrixV().col(5 + i);
      }

      return 0;
    }

    Err ComputeConstraintMatrix(const NullspaceBasis& basis,
                                std::vector<Polynomial33>& constraint) const
    {
      typedef hs::math::MultiPolynomial<Scalar, 3, 1> Polynomial31;
      typedef hs::math::MultiPolynomial<Scalar, 3, 2> Polynomial32;

      //由多项式矩阵的E矩阵，行主序排列
      std::vector<Polynomial31> E(9);
      for (size_t i = 0; i < 9; i++)
      {
        Polynomial31 e;
        e[0] = basis(i, 3);
        e[1] = basis(i, 0);
        e[2] = basis(i, 1);
        e[3] = basis(i, 2);
        E[i] = e;
      }

      //E矩阵充分必要条件约束
      std::vector<Polynomial32> EET(9);
      for (size_t i = 0; i < 3; i++)
      {
        for (size_t j = 0; j < 3; j++)
        {
          for (size_t k = 0; k < 3; k++)
          {
            EET[i * 3 + j] += E[i * 3 + k] * E[j * 3 + k];
          }
        }
      }

      Polynomial32 trace;
      for (size_t i = 0; i < 3; i++)
      {
        trace += EET[i * 3 + i];
      }

      //2EE^T-trace
      for (size_t i = 0; i < 3; i++)
      {
        for (size_t j = 0; j < 3; j++)
        {
          EET[i * 3 + j] *= 2;
          if (i == j)
          {
            EET[i * 3 + j] -= trace;
          }
        }
      }

      for (size_t i = 0; i < 3; i++)
      {
        for (size_t j = 0; j < 3; j++)
        {
          for (size_t k = 0; k < 3; k++)
          {
            constraint[i * 3 + j + 1] += EET[i * 3 + k] * E[k * 3 + j];
          }
        }
      }

      //E矩阵行列式为0的约束
      Polynomial33 determinant;
      determinant += E[0] * E[4] * E[8];
      determinant += E[1] * E[5] * E[6];
      determinant += E[2] * E[3] * E[7];
      determinant -= E[2] * E[4] * E[6];
      determinant -= E[1] * E[3] * E[8];
      determinant -= E[0] * E[5] * E[7];
      constraint[0] = determinant;

      return 0;
    }

    Err ComputeGroebnerBasis(const std::vector<Polynomial33>& constraint,
                             GroebnerBasis& groebner_basis) const
    {
      EIGEN_MATRIX(Scalar, 10, 20) A;

      //将多项式重新排列并放入A矩阵中
      //多项式顺序为
      //x3 x2y x2z xy2 xyz xz2 y3 y2z yz2 z3 x2 xy xz y2 yz z2 x y z 1
      for (size_t i = 0; i < 10; i++)
      {
        A(i, 0) = constraint[i][10];
        A(i, 1) = constraint[i][11];
        A(i, 2) = constraint[i][12];
        A(i, 3) = constraint[i][13];
        A(i, 4) = constraint[i][14];
        A(i, 5) = constraint[i][15];
        A(i, 6) = constraint[i][16];
        A(i, 7) = constraint[i][17];
        A(i, 8) = constraint[i][18];
        A(i, 9) = constraint[i][19];
        A(i, 10) = constraint[i][4];
        A(i, 11) = constraint[i][5];
        A(i, 12) = constraint[i][6];
        A(i, 13) = constraint[i][7];
        A(i, 14) = constraint[i][8];
        A(i, 15) = constraint[i][9];
        A(i, 16) = constraint[i][1];
        A(i, 17) = constraint[i][2];
        A(i, 18) = constraint[i][3];
        A(i, 19) = constraint[i][0];
      }

      //对A矩阵作高斯消去
      for (size_t i = 0; i < 10; i++)
      {
        Scalar leading = A(i, i);
        A.row(i) /= leading;
        for (size_t j = i + 1; j < 10; j++)
        {
          A.row(j) -= A.row(i) * A(j, i);
        }
      }
      for (size_t i = 9; i > 0; i--)
      {
        for (size_t j = 0; j < i; j++)
        {
          A.row(j) -= A.row(i) * A(j, i);
        }
      }

      groebner_basis = A.block(0, 10, 10, 10);

      return 0;
    }


    Err ComputeActionMatrix(const GroebnerBasis& groebner_basis,
                            ActionMatrix& action_matrix) const
    {
      action_matrix.setZero();
      action_matrix.row(0) = -groebner_basis.row(0);
      action_matrix.row(1) = -groebner_basis.row(1);
      action_matrix.row(2) = -groebner_basis.row(2);
      action_matrix.row(3) = -groebner_basis.row(4);
      action_matrix.row(4) = -groebner_basis.row(5);
      action_matrix.row(5) = -groebner_basis.row(7);
      action_matrix(6, 0) = 1;
      action_matrix(7, 1) = 1;
      action_matrix(8, 3) = 1;
      action_matrix(9, 6) = 1;

      return 0;
    }

    Err ComputeEMatrixHypothese(const ActionMatrix& action_matrix,
                                const NullspaceBasis& basis,
                                EMatrixHypotheses& ematrix_hypotheses) const
    {
      typedef Eigen::EigenSolver<ActionMatrix> ActionEigenSolver;
      typedef typename ActionEigenSolver::EigenvectorsType EigenVector;

      //计算action矩阵的实特征向量
      ActionEigenSolver eigen_solver(action_matrix);
      ematrix_hypotheses.clear();
      EigenVector eigen_vector = eigen_solver.eigenvectors();
      for (size_t i = 0; i < 10; i++)
      {
        if (eigen_solver.eigenvalues()[i].imag() == 0)
        {
          Scalar x = eigen_vector.col(i).real()[6];
          Scalar y = eigen_vector.col(i).real()[7];
          Scalar z = eigen_vector.col(i).real()[8];
          Scalar w = eigen_vector.col(i).real()[9];
          x /= w;
          y /= w;
          z /= w;

          EMatrix ematrix_hypothesis;
          for (size_t j = 0; j < 3; j++)
          {
            for (size_t k = 0; k < 3; k++)
            {
              ematrix_hypothesis(j, k) = x * basis.col(0)[j * 3 + k] +
                                         y * basis.col(1)[j * 3 + k] +
                                         z * basis.col(2)[j * 3 + k] +
                                             basis.col(3)[j * 3 + k];
            }
          }

          ematrix_hypothesis /= ematrix_hypothesis(2, 2);
          ematrix_hypotheses.push_back(ematrix_hypothesis);
        }
      }

      return 0;
    }
  };

  class EMatrixEvaluator
  {
  public:
    Err operator()(const HKeyPair& key_pair,
                   const EMatrix& ematrix,
                   Scalar& distance) const
    {
      HLine epiline_left = ematrix * key_pair.first;
      HLine epiline_right = ematrix.transpose() *
                            key_pair.second;

      Scalar error_left = std::abs(key_pair.second.dot(epiline_left)) /
                                    epiline_left.segment(0, 2).norm();
      Scalar error_right = std::abs(key_pair.first.dot(epiline_right)) /
                                    epiline_right.segment(0, 2).norm();
      distance = error_left + error_right;

      return 0;
    }
  };

  class HypothesesSelector
  {
  public:
    Err operator() (const HKeyPairContainer& key_pairs,
                    const EMatrixHypotheses& ematrix_hypotheses,
                    EMatrix& best_E) const
    {
      if (ematrix_hypotheses.empty())
      {
        return -1;
      }
      size_t number_of_hypotheses = ematrix_hypotheses.size();
      std::vector<Scalar> hypotheses_score(number_of_hypotheses, Scalar(0));
      EMatrixEvaluator evaluator;
      for (size_t i = 0; i < number_of_hypotheses; i++)
      {
        for (size_t j = 0; j < key_pairs.size(); j++)
        {
          Scalar distance;
          evaluator(key_pairs[j], ematrix_hypotheses[i], distance);
          hypotheses_score[i] += distance;
        }
      }

      Scalar min_score = std::numeric_limits<Scalar>::max();
      size_t best_id = 0;
      for (size_t i = 0; i < number_of_hypotheses; i++)
      {
        if (hypotheses_score[i] < min_score)
        {
          min_score = hypotheses_score[i];
          best_id = i;
        }
      }

      best_E = ematrix_hypotheses[best_id];

      return 0;
    }
  };

public:
  Err operator () (const HKeyPairContainer& key_pairs,
                   EMatrix& e_matrix) const
  {
    HypothesesGenerator hypotheses_generator;
    EMatrixHypotheses ematrix_hypotheses;
    if (hypotheses_generator(key_pairs, ematrix_hypotheses) != 0)
    {
      return -1;
    }

    HypothesesSelector hypotheses_selector;
    if (hypotheses_selector(key_pairs, ematrix_hypotheses, e_matrix) != 0)
    {
      return -1;
    }

    return 0;
  }
};

}
}
}

#endif
