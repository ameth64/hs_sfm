#ifndef _HS_SFM_UNIT_TEST_TEST_CAMERA_SHARED_NORMAL_EQUATION_SOLVER_HPP_
#define _HS_SFM_UNIT_TEST_TEST_CAMERA_SHARED_NORMAL_EQUATION_SOLVER_HPP_

#include <iostream>
#include <iomanip>
#include <limits>

#include "test_camera_shared_normal_equation_builder.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_augmentor.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_normal_equation_solver.hpp"

namespace test
{

template <typename _Scalar>
class TestCameraSharedNormalEquationSolver :
  public TestCameraSharedNormalEquationBuilder<_Scalar>
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef TestCameraSharedNormalEquationBuilder<Scalar> Base;

  typedef typename Base::VectorFunction VectorFunction;
  typedef typename Base::XVector XVector;
  typedef typename Base::YVector YVector;
  typedef typename Base::Index Index;
  typedef typename Base::NormalEquationBuilder NormalEquationBuilder;
  typedef typename Base::YCovarianceInverse YCovarianceInverse;
  typedef typename Base::NormalMatrix NormalMatrix;
  typedef typename Base::Gradient Gradient;
  typedef typename Base::Residuals Residuals;
  typedef typename Base::JacobianMatrixCalculator JacobianMatrixCalculator;
  typedef typename Base::JacobianMatrix JacobianMatrix;
  typedef typename Base::NoisedYGenerator NoisedYGenerator;
  typedef typename Base::KeyCovariance KeyCovariance;

  using Base::GenerateYCovarianceInverse;
  using Base::GenerateDenseYCovarianceInverse;
  using Base::CompareYCovarianceInverse;
  using Base::GenerateResiduals;
  using Base::CalculateJacobianMatrix;
  using Base::GenerateDenseJacobianMatrix;
  using Base::GenerateDenseNormalEquation;
  using Base::CalculateNormalEquation;
  using Base::CompareNormalEquation;

#ifdef IMAGE_TEST
  using Base::SaveMatrixImage;
#endif

  typedef hs::sfm::ba::CameraSharedNormalEquationSolver<Scalar>
    NormalEquationSolver;
  typedef typename NormalEquationSolver::DeltaXVector DeltaXVector;
  typedef hs::sfm::ba::CameraSharedAugmentor<Scalar> Augmentor;
  typedef typename Augmentor::AugmentedNormalMatrix AugmentedNormalMatrix;

private:
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, Eigen::Dynamic) MatrixXX;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) VectorX;

public:
  TestCameraSharedNormalEquationSolver(Scalar point_planar_stddev,
                                       Scalar point_height_stddev,
                                       Scalar image_rotation_stddev,
                                       Scalar image_position_stddev,
                                       Scalar camera_skew_stddev,
                                       Scalar camera_pixel_ratio_stddev)
    : Base(point_planar_stddev,
           point_height_stddev,
           image_rotation_stddev,
           image_position_stddev,
           camera_skew_stddev,
           camera_pixel_ratio_stddev)
  {

  }

  Err Test(const VectorFunction& vector_function,
           const XVector& x,
           const KeyCovariance& key_covariance) const
  {
    YCovarianceInverse y_covariance_inverse;
    if (GenerateYCovarianceInverse(vector_function,
                                   key_covariance,
                                   y_covariance_inverse) != 0)
    {
      return -1;
    }

    MatrixXX y_covariance_inverse_dense;
    if (GenerateDenseYCovarianceInverse(vector_function,
                                        key_covariance,
                                        y_covariance_inverse_dense) != 0)
    {
      return -1;
    }

    if (CompareYCovarianceInverse(y_covariance_inverse,
                                  y_covariance_inverse_dense) != 0)
    {
      return -1;
    }

    Residuals residuals;
    if (GenerateResiduals(vector_function,
                          x,
                          y_covariance_inverse,
                          residuals) != 0)
    {
      return -1;
    }

    JacobianMatrix jacobian_matrix;
    if (CalculateJacobianMatrix(vector_function, x, jacobian_matrix) != 0)
    {
      return -1;
    }

    MatrixXX jacobian_matrix_dense;
    if (GenerateDenseJacobianMatrix(jacobian_matrix,
                                    jacobian_matrix_dense) != 0)
    {
      return -1;
    }

    MatrixXX normal_matrix_dense;
    VectorX gradient_dense;
    if (GenerateDenseNormalEquation(jacobian_matrix_dense,
                                    residuals,
                                    y_covariance_inverse_dense,
                                    normal_matrix_dense,
                                    gradient_dense) != 0)
    {
      return -1;
    }

    NormalMatrix normal_matrix;
    Gradient gradient;
    if (CalculateNormalEquation(jacobian_matrix,
                                residuals,
                                y_covariance_inverse,
                                normal_matrix,
                                gradient) != 0)
    {
      return -1;
    }

    if (CompareNormalEquation(normal_matrix,
                              gradient,
                              normal_matrix_dense,
                              gradient_dense) != 0)
    {
      return -1;
    }

    Scalar mu = GenerateMu(normal_matrix);
    AugmentedNormalMatrix augmented_normal_matrix(normal_matrix, mu);

    MatrixXX augmented_normal_matrix_dense;
    if (GenerateAugmentedNormalMatrixDense(normal_matrix_dense, mu,
                                           augmented_normal_matrix_dense) != 0)
    {
      return -1;
    }

//#ifdef IMAGE_TEST
//    Index x_size = vector_function.GetXSize();
//    Index point_size = vector_function.GetPointParamsSize();
//    MatrixXX reduced_matrix_dense =
//      augmented_normal_matrix_dense.block(
//        point_size, point_size,
//        x_size - point_size, x_size - point_size) -
//      augmented_normal_matrix_dense.block(
//        point_size, 0,
//        x_size - point_size, point_size) *
//      augmented_normal_matrix_dense.block(0, 0,
//                                          point_size, point_size).inverse() *
//      augmented_normal_matrix_dense.block(
//        0, point_size,
//        point_size, x_size - point_size);
//    SaveMatrixImage(reduced_matrix_dense, "reduced_matrix_dense.png");
//    SaveMatrixImage(augmented_normal_matrix_dense,
//                    "augmented_normal_matrix_dense.png");
//    VectorX reduced_rhs_dense =
//      gradient_dense.segment(point_size, x_size - point_size) -
//      augmented_normal_matrix_dense.block(
//        point_size, 0,
//        x_size - point_size, point_size) *
//      augmented_normal_matrix_dense.block(0, 0,
//                                          point_size, point_size).inverse() *
//      gradient_dense.segment(0, point_size);
//
//    NormalEquationSolver solver;
//    typename NormalEquationSolver::ReducedMatrix reduced_matrix;
//    typename NormalEquationSolver::ReducedRHS reduced_rhs;
//    solver.ComputeReducedSystem(augmented_normal_matrix,
//                                gradient,
//                                reduced_matrix,
//                                reduced_rhs);
//
//    MatrixXX reduced_diff = MatrixXX::Zero(reduced_matrix_dense.rows(),
//                                           reduced_matrix_dense.cols());
//    std::cout.setf(std::ios::fixed);
//    std::cout<<std::setprecision(10);
//    for (Index i = 0; i < reduced_matrix_dense.rows(); i++)
//    {
//      for (Index j = 0; j < reduced_matrix_dense.cols(); j++)
//      {
//        Scalar value = reduced_matrix.coeff(i, j);
//        Scalar value_dense = reduced_matrix_dense.coeff(i, j);
//        Scalar error = std::abs(value - value_dense);
//        if (error > 1e-5)
//        {
//          std::cout<<"error:"<<error<<"\n";
//          std::cout<<"reduced_matrix_dense["<<i<<","<<j<<"]="<<value_dense<<"\n"
//                   <<"      reduced_matrix["<<i<<","<<j<<"]="<<value<<"\n";
//          reduced_diff(i, j) = error;
//        }
//      }
//    }
//    SaveMatrixImage(reduced_diff, "reduced_diff.png");
//
//    for (Index i = 0; i < reduced_rhs_dense.rows(); i++)
//    {
//      Scalar value = reduced_rhs[i];
//      Scalar value_dense = reduced_rhs_dense[i];
//      Scalar abs_error = std::abs(value - value_dense);
//      Scalar rel_error = abs_error;
//      if (value != Scalar(0)) rel_error = std::abs(rel_error / value);
//      if (abs_error > 1e-6 && rel_error > 1e-6)
//      {
//        std::cout<<"abs_error:"<<abs_error<<"\n"
//                 <<"rel_error:"<<rel_error<<"\n"
//                 <<"reduced_rhs_dense["<<i<<"]="<<value_dense<<"\n"
//                 <<"      reduced_rhs["<<i<<"]="<<value<<"\n";
//      }
//    }
//
//#endif

    DeltaXVector delta_x;
    if (SolveAugmentedNormalEquation(augmented_normal_matrix,
                                     gradient,
                                     delta_x) != 0)
    {
      return -1;
    }

    if (TestDeltaXVector(augmented_normal_matrix_dense, gradient_dense,
                         delta_x) != 0)
    {
      return -1;
    }

    return 0;
  }

  Scalar GenerateMu(const NormalMatrix& normal_matrix) const
  {
    Index x_size = normal_matrix.GetXSize();
    Scalar max_value = -std::numeric_limits<Scalar>::max();
    for (Index i = 0; i < x_size; i++)
    {
      max_value = std::max(max_value, normal_matrix.coeff(i, i));
    }
    return max_value;
  }

  Err GenerateAugmentedNormalMatrixDense(
    const MatrixXX& normal_matrix_dense,
    Scalar mu,
    MatrixXX& augmented_normal_matrix_dense) const
  {
    Index x_size = normal_matrix_dense.rows();
    augmented_normal_matrix_dense = normal_matrix_dense;
    for (Index i = 0; i < x_size; i++)
    {
      augmented_normal_matrix_dense(i, i) += mu;
    }
    return 0;
  }

  Err SolveAugmentedNormalEquation(
    const AugmentedNormalMatrix& augmented_normal_matrix,
    const Gradient& gradient,
    DeltaXVector& delta_x) const
  {
    NormalEquationSolver solver;
    return solver(augmented_normal_matrix, gradient, delta_x);
  }

  Err TestDeltaXVector(const MatrixXX& augmented_normal_matrix_dense,
                       const VectorX& gradient_dense,
                       const DeltaXVector& delta_x) const
  {
    Index x_size = delta_x.size();
    if (x_size != gradient_dense.rows() ||
        x_size != augmented_normal_matrix_dense.rows() ||
        x_size != augmented_normal_matrix_dense.cols())
    {
      return -1;
    }
    VectorX gradient = augmented_normal_matrix_dense * delta_x;
    const Scalar threshold = 1e-5;
    Err result = 0;
    for (Index i = 0; i < x_size; i++)
    {
      Scalar value = gradient[i];
      Scalar value_dense = gradient_dense[i];
      Scalar error_abs = std::abs(value - value_dense);
      Scalar error_rel = error_abs;
      if (value_dense != 0) error_rel /= std::abs(value_dense);
      if (error_abs > threshold && error_rel > threshold)
      {
        std::cout<<"error_abs:"<<error_abs<<"\n";
        std::cout<<"error_rel:"<<error_rel<<"\n";
        std::cout<<"      gradient["<<i<<"]="<<value<<"\n";
        std::cout<<"gradient_dense["<<i<<"]="<<value_dense<<"\n";
        result = -1;
      }
    }

    return result;
  }
};

}

#endif
