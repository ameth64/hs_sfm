#ifndef _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_DLT_HPP_
#define _HS_SFM_TRIANGULATE_MULTIPLE_VIEW_DLT_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace triangulate
{

template <typename _Scalar>
class MultipleViewDLT
{
public:
  typedef _Scalar Scalar;
  typedef EIGEN_MAT(Scalar, 3, 4) PMatrix;
  typedef EIGEN_VECTOR(PMatrix) PMatrixContainer;
  typedef EIGEN_VEC(Scalar, 3) HomogeneousKey;
  typedef EIGEN_VECTOR(HomogeneousKey) HomogeneousKeyContainer;
  typedef EIGEN_VEC(Scalar, 3) Point;
  typedef EIGEN_VEC(Scalar, 4) HomogeneousPoint;
  typedef EIGEN_MAT(Scalar, Eigen::Dynamic, 4) AMatrix;
  typedef typename PMatrix::Index Index;

  typedef int Err;

  Err operator()(const PMatrixContainer& p_matrices,
                 const HomogeneousKeyContainer& homogeneous_keys,
                 Point& point) const
  {
    Index number_of_views = Index(p_matrices.size());
    if (number_of_views != Index(homogeneous_keys.size()))
    {
      return -1;
    }
    AMatrix A;
    A.resize(2 * number_of_views, 4);
    for (Index i = 0; i < number_of_views; i++)
    {
      A.row(i * 2 + 0) = homogeneous_keys[i][0] * p_matrices[i].row(2) -
                         p_matrices[i].row(0);
      A.row(i * 2 + 1) = homogeneous_keys[i][1] * p_matrices[i].row(2) -
                         p_matrices[i].row(1);
    }

    Eigen::JacobiSVD<AMatrix> svd(A, Eigen::ComputeFullV);
    HomogeneousPoint homogeneous_point = svd.matrixV().col(3);
    point = homogeneous_point.segment(0, 3) / homogeneous_point[3];
    return 0;
  }
};

}
}
}

#endif
