#ifndef _HS_SFM_HOMOGRAPHY_HOMOGRAPHY2D_ANALYTICAL_JACOBIAN_MATRIX_CALCULATOR_HPP_
#define _HS_SFM_HOMOGRAPHY_HOMOGRAPHY2D_ANALYTICAL_JACOBIAN_MATRIX_CALCULATOR_HPP_

#include "hs_sfm/homography/homography2d_vector_function.hpp"
#include "hs_sfm/homography/homography2d_jacobian_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace homography
{

template <typename _VectorFunction>
class Homography2DAnalyticalJacobianMatrixCalculator;

template <typename _Scalar>
class Homography2DAnalyticalJacobianMatrixCalculator<
        Homography2DVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef Homography2DVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef Homography2DJacobianMatrix<Scalar> JacobianMatrix;

private:
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::HomogeneousMatrix HomogeneousMatrix;
  typedef typename VectorFunction::Key Key;
  typedef typename VectorFunction::HomogeneousKey HomogeneousKey;
  typedef typename JacobianMatrix::KeyBlock KeyBlock;
  typedef typename JacobianMatrix::KeyBlockContainer KeyBlockContainer;
  typedef typename JacobianMatrix::HBlock HBlock;
  typedef typename JacobianMatrix::HBlockContainer HBlockContainer;
  
public:
  Scalar m_a;
  Err operator() (const VectorFunction& vector_function,
                  const XVector& x,
                  JacobianMatrix& jacobian_matrix) const
  {
    Index x_size = vector_function.GetXSize();
    if (x.rows() != x_size)
    {
      return -1;
    }
    Index y_size = vector_function.GetYSize();
    Index number_of_keys = vector_function.number_of_keys();
    HomogeneousMatrix H;
    vector_function.GetHomographyMatrix(x, H);
    KeyBlockContainer& key_blocks = jacobian_matrix.key_blocks();
    key_blocks.resize(number_of_keys);
    HBlockContainer& h_blocks = jacobian_matrix.h_blocks();
    h_blocks.resize(number_of_keys);
    for (Index i = 0; i < number_of_keys; i++)
    {
       Key key;
       vector_function.GetXKey(x, i, key);
       Scalar w = H.template block<1, 2>(2, 0) * key + H(2, 2);
       Key transformed_key = (H.template block<2, 2>(0, 0) * key +
                              H.template block<2, 1>(0, 2)) / w;
       key_blocks[i].row(0) = (H.template block<1, 2>(0, 0) -
                               transformed_key[0] *
                               H.template block<1, 2>(2, 0)) / w;
       key_blocks[i].row(1) = (H.template block<1, 2>(1, 0) -
                               transformed_key[1] *
                               H.template block<1, 2>(2, 0)) / w;
       h_blocks[i].template block<1, 2>(0, 0) = key.transpose() / w;
       h_blocks[i](0, 2) = Scalar(1) / w;
       h_blocks[i].template block<1, 3>(0, 3).setZero();
       h_blocks[i].template block<1, 2>(0, 6) =
         -transformed_key[0] * key.transpose() / w;
       h_blocks[i](0, 8) = -transformed_key[0] / w;
       h_blocks[i].template block<1, 3>(1, 0).setZero();
       h_blocks[i].template block<1, 2>(1, 3) = key.transpose() / w;
       h_blocks[i](1, 5) = Scalar(1) / w;
       h_blocks[i].template block<1, 2>(1, 6) =
         -transformed_key[1] * key.transpose() / w;
       h_blocks[i](1, 8) = -transformed_key[1] / w;
    }
    
    return 0;
  }
};

}
}
}

#endif