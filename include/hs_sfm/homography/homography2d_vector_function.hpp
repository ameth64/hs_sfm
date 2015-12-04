#ifndef _HS_SFM_HOMOGRAPHY_HOMOGRAPHY2D_VECTOR_FUNCTION_HPP_
#define _HS_SFM_HOMOGRAPHY_HOMOGRAPHY2D_VECTOR_FUNCTION_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace homography
{

template <typename _Scalar>
class Homography2DVectorFunction
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) XVector;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) YVector;
  typedef typename XVector::Index Index;
  typedef EIGEN_MATRIX(Scalar, 3, 3) HomogeneousMatrix;
  typedef EIGEN_VECTOR(Scalar, 2) Key;
  typedef EIGEN_VECTOR(Scalar, 3) HomogeneousKey;

private:
  enum
  {
    Major = Eigen::EIGEN_DEFAULT_MATRIX_STORAGE_ORDER_OPTION
  };

public:
  Homography2DVectorFunction() : number_of_keys_(0) {}
  Homography2DVectorFunction(Index number_of_keys)
    : number_of_keys_(number_of_keys) {}

  Err operator() (const XVector& x, YVector& y) const
  {
    Index x_size = GetXSize();
    Index y_size = GetYSize();
    if (x.rows() != x_size)
    {
      return -1;
    }
    y.resize(y_size);
    HomogeneousMatrix H;
    GetHomographyMatrix(x, H);
    for (Index i = 0; i < number_of_keys_; i++)
    {
      HomogeneousKey transformed_hkey =
        H.template block<3, 2>(0, 0) *
        x.template segment<2>(i * 2) +
        H.col(2);
      transformed_hkey /= transformed_hkey[2];
      y.template segment<2>(i * 4) = x.template segment<2>(i * 2);
      y.template segment<2>(i * 4 + 2) =
        transformed_hkey.template segment<2>(0);
    }

    return 0;
  }

  inline Index number_of_keys() const {return number_of_keys_;}

  inline void set_number_of_keys(Index number_of_keys)
  {
    number_of_keys_ = number_of_keys;
  }

  inline Index GetXSize() const
  {
    return number_of_keys_ * 2 + 9;
  }

  inline Index GetYSize() const
  {
    return number_of_keys_ * 4;
  }

  Err GetHomographyMatrix(const XVector& x,
                          HomogeneousMatrix& homogeneous_matrix) const
  {
    Index x_size = GetXSize();
    if (x.rows() != x_size)
    {
      return -1;
    }

    Index start_id = x_size - 9;
    //if (Major == Eigen::ColMajor)
    //{
    //  homogeneous_matrix << x[start_id + 0], x[start_id + 3], x[start_id + 6],
    //                        x[start_id + 1], x[start_id + 4], x[start_id + 7],
    //                        x[start_id + 2], x[start_id + 5], x[start_id + 8];
    //}
    //else
    {
      homogeneous_matrix << x[start_id + 0], x[start_id + 1], x[start_id + 2],
                            x[start_id + 3], x[start_id + 4], x[start_id + 5],
                            x[start_id + 6], x[start_id + 7], x[start_id + 8];
    }

    return 0;
  }

  Err SetHomographyMatrix(const HomogeneousMatrix& homo_geneous_matrix,
                          XVector& x) const
  {
    Index x_size = GetXSize();
    if (x.rows() != x_size)
    {
      return -1;
    }

    Index start_id = x_size - 9;
    //if (Major == Eigen::ColMajor)
    //{
    //  x[start_id + 0] = homo_geneous_matrix(0, 0);
    //  x[start_id + 1] = homo_geneous_matrix(1, 0);
    //  x[start_id + 2] = homo_geneous_matrix(2, 0);
    //  x[start_id + 3] = homo_geneous_matrix(0, 1);
    //  x[start_id + 4] = homo_geneous_matrix(1, 1);
    //  x[start_id + 5] = homo_geneous_matrix(2, 1);
    //  x[start_id + 6] = homo_geneous_matrix(0, 2);
    //  x[start_id + 7] = homo_geneous_matrix(1, 2);
    //  x[start_id + 8] = homo_geneous_matrix(2, 2);
    //}
    //else
    {
      x[start_id + 0] = homo_geneous_matrix(0, 0);
      x[start_id + 1] = homo_geneous_matrix(0, 1);
      x[start_id + 2] = homo_geneous_matrix(0, 2);
      x[start_id + 3] = homo_geneous_matrix(1, 0);
      x[start_id + 4] = homo_geneous_matrix(1, 1);
      x[start_id + 5] = homo_geneous_matrix(1, 2);
      x[start_id + 6] = homo_geneous_matrix(2, 0);
      x[start_id + 7] = homo_geneous_matrix(2, 1);
      x[start_id + 8] = homo_geneous_matrix(2, 2);
    }

    return 0;
  }

  Err GetXKey(const XVector& x, Index i, Key& key) const
  {
    Index x_size = GetXSize();
    if (x.rows() != x_size)
    {
      return -1;
    }

    key = x.template segment<2>(i * 2);

    return 0;
  }

  Err GetYKey(const YVector& y, Index i, Key& key) const
  {
    Index y_size = GetYSize();
    if (y.rows() != y_size)
    {
      return -1;
    }

    key = y.template segment<2>(i * 4);

    return 0;
  }

  Err GetYTransformedKey(const YVector& y, Index i, Key& transformed_key) const
  {
    Index y_size = GetYSize();
    if (y.rows() != y_size)
    {
      return -1;
    }

    transformed_key = y.template segment<2>(i * 4 + 2);

    return 0;
  }

private:
  Index number_of_keys_;
};

}
}
}

#endif
