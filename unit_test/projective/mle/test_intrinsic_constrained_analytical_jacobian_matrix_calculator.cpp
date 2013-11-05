#include <iostream>

#include <gtest/gtest.h>

#include "hs_math/linear_algebra/latraits/vector_eigen.hpp"
#include "hs_math/linear_algebra/lafunc/arithmetic_eigen.hpp"
#include "hs_math/fdjac/forward_finite_difference_jacobian_matrix_calculator.hpp"

#include "hs_sfm/projective/mle/intrinsic_constrained_synthetic_generator.hpp"
#include "hs_sfm/projective/mle/intrinsic_constrained_vector_function.hpp"

#include "hs_sfm/projective/mle/intrinsic_constrained_analytical_jacobian_matrix_calculator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestIntrinsicConstrainedAnalyticalJacobianMatrixCalculator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

private:
  typedef hs::sfm::projective::IntrinsicConstrainedSyntheticGenerator<
            Scalar, ImageDimension> SyntheticGenerator;
  typedef hs::sfm::projective::IntrinsicConstrainedVectorFunction<Scalar>
          VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef
    hs::sfm::projective::IntrinsicConstrainedAnalyticalJacobianMatrixCalculator<
      VectorFunction> AnalyticalJacobianMatrixCalculator;
  typedef typename AnalyticalJacobianMatrixCalculator::JacobianMatrix
          AnalyticalJacobianMatrix;
  typedef
    hs::math::fdjac::ForwardFiniteDifferenceJacobianMatrixCalculator<
      VectorFunction> FFDJacobianMatrixCalculator;
  typedef typename FFDJacobianMatrixCalculator::JacobianMatrix
          FFDJacobianMatrix;

public:
  TestIntrinsicConstrainedAnalyticalJacobianMatrixCalculator(
    Scalar focal_length_in_metre,
    Scalar ground_resolution,
    ImageDimension image_width,
    ImageDimension image_height,
    Scalar pixel_size,
    size_t number_of_points,
    Scalar scene_max_height,
    Scalar camera_height_stddev,
    Scalar camera_planar_stddev,
    Scalar camera_rot_stddev,
    Scalar north_west_angle,
    Scalar skew,
    Scalar principal_point_x,
    Scalar principal_point_y,
    Scalar pixel_ratio)
    : synthetic_generator_(focal_length_in_metre,
                           ground_resolution,
                           image_width,
                           image_height,
                           pixel_size,
                           number_of_points,
                           scene_max_height,
                           camera_height_stddev,
                           camera_planar_stddev,
                           camera_rot_stddev,
                           north_west_angle,
                           skew,
                           principal_point_x,
                           principal_point_y,
                           pixel_ratio) {}

  Err Test()
  {
    VectorFunction vector_function;
    XVector x;
    YVector y;
    if (synthetic_generator_(vector_function, x, y) != 0)
    {
      return -1;
    }

    AnalyticalJacobianMatrixCalculator analytical_calculator;
    AnalyticalJacobianMatrix analytical_jacobian_matrix;
    if (analytical_calculator(vector_function, x,
                              analytical_jacobian_matrix) != 0)
    {
      return -1;
    }

    FFDJacobianMatrixCalculator ffd_calculator(Scalar(1e-6), Scalar(1e-8));
    FFDJacobianMatrix ffd_jacobian_matrix;
    if (ffd_calculator(vector_function, x,
                       ffd_jacobian_matrix) != 0)
    {
      return -1;
    }

    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    Err result = 0;
    for (Index i = 0; i < y_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        Scalar ffd_value = ffd_jacobian_matrix.coeff(i, j);
        Scalar analytical_value = analytical_jacobian_matrix.coeff(i, j);
        Scalar absolute_error = abs(analytical_value - ffd_value);
        Scalar relative_error = absolute_error;
        if (ffd_value != Scalar(0)) relative_error /= ffd_value;
        if (relative_error > Scalar(5e-4) && absolute_error > Scalar(1e-4))
        {
          std::cout<<"analytical_jacobian_matrix[i, j]="<<analytical_value<<".\n";
          std::cout<<"ffd_jacobian_matrix[i, j]="<<ffd_value<<".\n";
          result = -1;
        }
      }
    }

    return result;
  }

private:
  SyntheticGenerator synthetic_generator_;
};

TEST(TestIntrinsicConstrainedAnalyticalJacobianMatrixCalculator, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestIntrinsicConstrainedAnalyticalJacobianMatrixCalculator<
            Scalar, ImageDimension> Test;

  Scalar focal_length_in_metre = 0.019;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 10;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rot_stddev = 1;
  Scalar north_west_angle = 60;
  Scalar skew = 0.0001;
  Scalar principal_point_x = 20.0;
  Scalar principal_point_y = 30.0;
  Scalar pixel_ratio = 1.0001;

  Test test(focal_length_in_metre,
            ground_resolution,
            image_width,
            image_height,
            pixel_size,
            number_of_points,
            scene_max_height,
            camera_height_stddev,
            camera_planar_stddev,
            camera_rot_stddev,
            north_west_angle,
            skew,
            principal_point_x,
            principal_point_y,
            pixel_ratio);

  ASSERT_EQ(0, test.Test());
}

}