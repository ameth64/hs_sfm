#include <iostream>
#include <ctime>

#include <gtest/gtest.h>

#include "hs_math/linear_algebra/latraits/vector_eigen.hpp"
#include "hs_math/linear_algebra/lafunc/arithmetic_eigen.hpp"
#include "hs_math/fdjac/forward_finite_difference_sparse_jacobian_matrix_calculator.hpp"

#include "hs_sfm/homography/homography2d_vector_function.hpp"
#include "hs_sfm/homography/homography2d_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/homography/homography2d_synthetic_data_generator.hpp"

namespace
{

template <typename _Scalar>
class TestHomography2DAnalyticalJacobianMatrixCalculator
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::homography::Homography2DVectorFunction<Scalar>
          VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;

  typedef
    hs::math::fdjac::ForwardFiniteDifferenceSparseJacobianMatrixCalculator<
      VectorFunction>
    FFDJacobianMatrixCalculator;
  typedef typename FFDJacobianMatrixCalculator::JacobianMatrix
                   FFDJacobianMatrix;

  typedef hs::sfm::homography::Homography2DAnalyticalJacobianMatrixCalculator<
            VectorFunction>
          AnalyticalJacobianMatrixCalculator;
  typedef typename AnalyticalJacobianMatrixCalculator::JacobianMatrix
                   AnalyticalJacobianMatrix;

public:
  static Err Test(const VectorFunction& vector_function,
                  const XVector& x)
  {
    FFDJacobianMatrixCalculator ffd_calculator(Scalar(1e-6),
                                               Scalar(1e-8),
                                               Scalar(1e-12));
    AnalyticalJacobianMatrixCalculator analytical_calculator;

    FFDJacobianMatrix ffd_jacobian_matrix;
    AnalyticalJacobianMatrix analytical_jacobian_matrix;

    time_t ffd_begin = clock();
    if (0 != ffd_calculator(vector_function, x, ffd_jacobian_matrix)) return -1;
    time_t ffd_end = clock();
    std::cout<<"FFD took "
              <<Scalar(ffd_end - ffd_begin) /
                Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";

    time_t analytical_begin = clock();
    if (0 != analytical_calculator(vector_function, x,
                                   analytical_jacobian_matrix)) return -1;
    time_t analytical_end = clock();
    std::cout<<"analytical took "
             <<Scalar(analytical_end - analytical_begin) /
               Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";

    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();
    Err result = 0;
    for (Index i = 0; i < y_size; i++)
    {
      for (Index j = 0; j < x_size; j++)
      {
        Scalar ffd_value = ffd_jacobian_matrix.coeff(i, j);
        Scalar analytical_value = analytical_jacobian_matrix.coeff(i, j);
        Scalar abs_error = std::abs(analytical_value - ffd_value);
        Scalar rel_error = (analytical_value == 0 ? abs_error :
                                                    abs_error / analytical_value);
        if (rel_error > Scalar(1e-4))
        {
          std::cout<<"analytical_jacobian_matrix[i, j]="
                   <<analytical_value<<".\n";
          std::cout<<"ffd_jacobian_matrix[i, j]="
                   <<ffd_value<<".\n";
          result = -1;
        }
      }
    }

    return 0;
  }
};

TEST(TestHomography2DAnalyticalJacobianMatrixCalculator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef hs::sfm::homography::Homography2DSyntheticDataGenerator<
            Scalar, ImageDimension> DataGenerator;
  typedef hs::sfm::homography::Homography2DVectorFunction<Scalar>
          VectorFunction;
  typedef VectorFunction::XVector XVector;
  typedef VectorFunction::YVector YVector;
  typedef TestHomography2DAnalyticalJacobianMatrixCalculator<Scalar> Test;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 1;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 1;
  Scalar north_west_angle = 60;
  size_t number_of_keys = 15;
  Scalar points_height = 25;

  DataGenerator data_generator(focal_length_in_metre,
                               number_of_strips,
                               ground_resolution,
                               image_width,
                               image_height,
                               pixel_size,
                               lateral_overlap_ratio,
                               longitudinal_overlap_ratio,
                               scene_max_height,
                               camera_height_stddev,
                               camera_planar_stddev,
                               camera_rotation_stddev,
                               north_west_angle,
                               number_of_keys,
                               points_height);

  VectorFunction vector_function;
  XVector x;
  YVector y;
  ASSERT_EQ(0, data_generator(vector_function, x, y));
  ASSERT_EQ(0, Test::Test(vector_function, x));
}

}