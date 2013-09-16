#include <iostream>
#include <ctime>

#include <gtest/gtest.h>

#include "hs_math/linear_algebra/latraits/vector_eigen.hpp"
#include "hs_math/linear_algebra/lafunc/arithmetic_eigen.hpp"
#include "hs_math/fdjac/forward_finite_difference_sparse_jacobian_matrix_calculator.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_forward_finite_difference_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_synthetic_data_generator.hpp"

namespace
{

template <typename _Scalar>
class TestBANaiveJacobianMatrixCalculator
{
public:
  typedef _Scalar Scalar;

  typedef int Err;

  typedef hs::sfm::ba::BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;

  typedef
    hs::math::fdjac::ForwardFiniteDifferenceSparseJacobianMatrixCalculator<
      VectorFunction>
    FFDJacobianMatrixCalculator;
  typedef typename FFDJacobianMatrixCalculator::JacobianMatrix
                   FFDJacobianMatrix;

  typedef hs::sfm::ba::BANaiveAnalyticalJacobianMatrixCalculator<
            VectorFunction>
          BAAnalyticalJacobianMatrixCalculator;
  typedef typename BAAnalyticalJacobianMatrixCalculator::JacobianMatrix
                   BAAnaliticalJacobianMatrix;

  typedef
    hs::sfm::ba::BANaiveForwardFiniteDifferenceJacobianMatrixCalculator<
      VectorFunction>
    BAFFDJacobianMatrixCalculator;
  typedef typename BAFFDJacobianMatrixCalculator::JacobianMatrix
                   BAFFDJacobianMatrix;

  static Err Test(const VectorFunction& vector_function,
                  const XVector& x, bool with_ffd)
  {
    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();

    FFDJacobianMatrixCalculator ffd_calculator(
      Scalar(1e-6), Scalar(1e-8), Scalar(1e-10));
    BAAnalyticalJacobianMatrixCalculator ba_analytical_calculator;
    BAFFDJacobianMatrixCalculator ba_ffd_calculator(Scalar(1e-6), Scalar(1e-8));

    FFDJacobianMatrix ffd_jacobian_matrix;
    BAAnaliticalJacobianMatrix ba_analytical_jacobian_matrix;
    BAFFDJacobianMatrix ba_ffd_jacobian_matrix;

    if (with_ffd)
    {
      time_t ffd_begin = clock();
      if (0 != ffd_calculator(vector_function, x, ffd_jacobian_matrix)) return -1;
      time_t ffd_end = clock();
      std::cout<<"FFD took "
               <<Scalar(ffd_end - ffd_begin) /
                 Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";
    }
    time_t ba_analytical_begin = clock();
    if (0 != ba_analytical_calculator(vector_function, x, ba_analytical_jacobian_matrix)) return -1;
    time_t ba_analytical_end = clock();
    std::cout<<"BA Analytic took "
             <<Scalar(ba_analytical_end - ba_analytical_begin) /
               Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";
    time_t ba_ffd_begin = clock();
    if (0 != ba_ffd_calculator(vector_function, x, ba_ffd_jacobian_matrix)) return -1;
    time_t ba_ffd_end = clock();
    std::cout<<"BA FFD took "
             <<Scalar(ba_ffd_end - ba_ffd_begin) /
               Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";

    if (ba_analytical_jacobian_matrix.number_of_cameras() !=
          ba_ffd_jacobian_matrix.number_of_cameras() ||
        ba_analytical_jacobian_matrix.number_of_points() !=
          ba_ffd_jacobian_matrix.number_of_points() ||
        ba_analytical_jacobian_matrix.camera_derivatives().size() !=
          ba_ffd_jacobian_matrix.camera_derivatives().size() ||
        ba_analytical_jacobian_matrix.point_derivatives().size() !=
          ba_ffd_jacobian_matrix.point_derivatives().size())
    {
      std::cout<<"Number of cameras or number of points not match.\n";
      return -1;
    }

    Err result = 0;
    const Scalar threshold = Scalar(1e-5);
    for (Index i = 0;
         i < Index(ba_ffd_jacobian_matrix.camera_derivatives().size());
         i++)
    {
      const typename BAAnaliticalJacobianMatrix::CameraDerivativeBlock&
        ba_analytical_camera_block =
          ba_analytical_jacobian_matrix.camera_derivatives()[i];
      const typename BAFFDJacobianMatrix::CameraDerivativeBlock&
        ba_ffd_camera_block =
          ba_ffd_jacobian_matrix.camera_derivatives()[i];

      if (ba_analytical_camera_block.derivative_block.rows() != 
          ba_ffd_camera_block.derivative_block.rows() ||
          ba_analytical_camera_block.derivative_block.cols() != 
          ba_ffd_camera_block.derivative_block.cols())
      {
        std::cout<<"Size of camera derivative matrix not match.\n";
        result = -1;
        break;
      }

      for (Index j = 0; j < ba_ffd_camera_block.derivative_block.rows(); j++)
      {
        for (Index k = 0; k < ba_ffd_camera_block.derivative_block.cols(); k++)
        {
          Scalar ba_ffd_value =
            ba_ffd_camera_block.derivative_block(j, k);
          Scalar ba_analytical_value =
            ba_analytical_camera_block.derivative_block(j, k);
          Scalar abs_err = std::abs(ba_ffd_value - ba_analytical_value);
          if (abs_err > threshold)
          {
            std::cout<<"Deference between the camera derivative matrix value is too big\n";
            std::cout<<"Camera derivative matrix Id:"<<i<<" .\n";
            std::cout<<"ba_ffd_camera_block["<<j<<", "<<k<<"] = "
                     <<ba_ffd_value<<".\n";
            std::cout<<"ba_analytical_camera_block["<<j<<", "<<k<<"] = "
                     <<ba_analytical_value<<".\n";
            result = -1;
          }
        }
      }
    }

    for (Index i = 0;
         i < Index(ba_ffd_jacobian_matrix.point_derivatives().size());
         i++)
    {
      const typename BAAnaliticalJacobianMatrix::PointDerivativeBlock&
        ba_analytical_point_block =
          ba_analytical_jacobian_matrix.point_derivatives()[i];
      const typename BAFFDJacobianMatrix::PointDerivativeBlock&
        ba_ffd_point_block =
          ba_ffd_jacobian_matrix.point_derivatives()[i];

      if (ba_analytical_point_block.derivative_block.rows() != 
          ba_ffd_point_block.derivative_block.rows() ||
          ba_analytical_point_block.derivative_block.cols() != 
          ba_ffd_point_block.derivative_block.cols())
      {
        std::cout<<"Size of point derivative matrix not match.\n";
        result = -1;
        break;
      }

      for (Index j = 0; j < ba_ffd_point_block.derivative_block.rows(); j++)
      {
        for (Index k = 0; k < ba_ffd_point_block.derivative_block.cols(); k++)
        {
          Scalar ba_ffd_value =
            ba_ffd_point_block.derivative_block(j, k);
          Scalar ba_analytical_value =
            ba_analytical_point_block.derivative_block(j, k);
          Scalar abs_err = std::abs(ba_ffd_value - ba_analytical_value);
          if (abs_err > threshold)
          {
            std::cout<<"Deference between the point derivative matrix value is too big\n";
            std::cout<<"Point derivative matrix Id:"<<i<<" .\n";
            std::cout<<"ba_ffd_point_block["<<j<<", "<<k<<"] = "
                     <<ba_ffd_value<<".\n";
            std::cout<<"ba_analytical_point_block["<<j<<", "<<k<<"] = "
                     <<ba_analytical_value<<".\n";
            result = -1;
          }
        }
      }
    }

    if (with_ffd)
    {
      for (Index i = 0; i < y_size; i++)
      {
        for (Index j = 0; j < x_size; j++)
        {
          Scalar ffd_value = ffd_jacobian_matrix.coeff(i, j);
          Scalar ba_ffd_value = ba_ffd_jacobian_matrix.coeff(i, j);
          Scalar abs_error = abs(ba_ffd_value - ffd_value);
          if (abs_error > Scalar(1e-8))
          {
            std::cout<<"ba_ffd_jacobian_matrix[i, j]="<<ba_ffd_value<<".\n";
            std::cout<<"ffd_jacobian_matrix[i, j]="<<ffd_value<<".\n";
            return -1;
          }
        }
      }
    }

    return result;
  }
};

TEST(TestBANaiveJacobianMatrixCalculator, SmallDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef hs::sfm::ba::BANaiveSyntheticDataGenerator<Scalar, ImageDimension>
          DataGenerator;
  typedef TestBANaiveJacobianMatrixCalculator<Scalar> Test;
  typedef Test::VectorFunction VectorFunction;
  typedef VectorFunction::XVector XVector;
  typedef VectorFunction::YVector YVector;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 3;
  size_t number_of_cameras_in_strip = 3;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 30;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 1;
  Scalar north_west_angle = 60;

  VectorFunction vector_function;
  XVector x;
  YVector y;

  DataGenerator data_generator(focal_length_in_metre,
                               number_of_strips,
                               number_of_cameras_in_strip,
                               ground_resolution,
                               image_width,
                               image_height,
                               pixel_size,
                               number_of_points,
                               lateral_overlap_ratio,
                               longitudinal_overlap_ratio,
                               scene_max_height,
                               camera_height_stddev,
                               camera_planar_stddev,
                               camera_rotation_stddev,
                               north_west_angle);
  ASSERT_EQ(0, data_generator(vector_function, x, y));

  ASSERT_EQ(0, Test::Test(vector_function, x, true));
}

TEST(TestBANaiveJacobianMatrixCalculator, BigDataTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef hs::sfm::ba::BANaiveSyntheticDataGenerator<Scalar, ImageDimension> DataGenerator;
  typedef TestBANaiveJacobianMatrixCalculator<Scalar> Test;
  typedef Test::VectorFunction VectorFunction;
  typedef VectorFunction::XVector XVector;
  typedef VectorFunction::YVector YVector;

  Scalar focal_length_in_metre = 0.019;
  size_t number_of_strips = 10;
  size_t number_of_cameras_in_strip = 10;
  Scalar ground_resolution = 0.1;
  ImageDimension image_width = 6000;
  ImageDimension image_height = 4000;
  Scalar pixel_size = 0.0000039;
  size_t number_of_points = 2000;
  Scalar lateral_overlap_ratio = 0.6;
  Scalar longitudinal_overlap_ratio = 0.8;
  Scalar scene_max_height = 50;
  Scalar camera_height_stddev = 5;
  Scalar camera_planar_stddev = 5;
  Scalar camera_rotation_stddev = 1;
  Scalar north_west_angle = 60;

  VectorFunction vector_function;
  XVector x;
  YVector y;

  DataGenerator data_generator(focal_length_in_metre,
                               number_of_strips,
                               number_of_cameras_in_strip,
                               ground_resolution,
                               image_width,
                               image_height,
                               pixel_size,
                               number_of_points,
                               lateral_overlap_ratio,
                               longitudinal_overlap_ratio,
                               scene_max_height,
                               camera_height_stddev,
                               camera_planar_stddev,
                               camera_rotation_stddev,
                               north_west_angle);
  ASSERT_EQ(0, data_generator(vector_function, x, y));

  ASSERT_EQ(0, Test::Test(vector_function, x, false));

}

}
