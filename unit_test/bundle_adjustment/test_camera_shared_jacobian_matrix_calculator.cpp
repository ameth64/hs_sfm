#include <iostream>
#include <iomanip>
#include <ctime>

#include "gtest/gtest.h"

#include "hs_math/linear_algebra/latraits/vector_eigen.hpp"
#include "hs_math/linear_algebra/lafunc/arithmetic_eigen.hpp"
#include "hs_math/fdjac/forward_finite_difference_sparse_jacobian_matrix_calculator.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_synthetic_data_generator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_analytical_jacobian_matrix_calculator.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_forward_finite_difference_jacobian_matrix_calculator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestCameraSharedJacobianMatrixCalculator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

  typedef hs::sfm::ba::CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;

  typedef
    hs::math::fdjac::ForwardFiniteDifferenceSparseJacobianMatrixCalculator<
      VectorFunction>
    FFDJacobianMatrixCalculator;
  typedef typename FFDJacobianMatrixCalculator::JacobianMatrix
                   FFDJacobianMatrix;

  typedef
    hs::sfm::ba::CameraSharedForwardFiniteDefferenceJacobianMatrixCalculator<
      VectorFunction>
    BAFFDJacobianMatrixCalculator;
  typedef typename BAFFDJacobianMatrixCalculator::JacobianMatrix
                   BAFFDJacobianMatrix;

  typedef hs::sfm::ba::CameraSharedAnalyticalJacobianMatrixCalculator<
            VectorFunction>
          BAAnalitycalJacobianMatrixCalculator;
  typedef typename BAAnalitycalJacobianMatrixCalculator::JacobianMatrix
                   BAAnalitycalJacobianMatrix;

public:
  Err Test(const VectorFunction& vector_function,
           const XVector& x, bool with_ffd) const
  {
    Index x_size = vector_function.GetXSize();
    Index y_size = vector_function.GetYSize();

    Scalar delta = Scalar(1e-6);
    Scalar min_delta = Scalar(1e-8);
    Scalar zero_threshold = Scalar(1e-10);
    FFDJacobianMatrixCalculator ffd_calculator(delta,
                                               min_delta,
                                               zero_threshold);
    BAFFDJacobianMatrixCalculator ba_ffd_calculator(delta, min_delta);
    BAAnalitycalJacobianMatrixCalculator ba_analytical_calculator;

    FFDJacobianMatrix ffd_jacobian_matrix;
    BAFFDJacobianMatrix ba_ffd_jacobian_matrix;
    BAAnalitycalJacobianMatrix ba_analytical_jacobian_matrix;

    std::cout.setf(std::ios::fixed);
    std::cout<<std::setprecision(10);

    if (with_ffd)
    {
      time_t ffd_begin = clock();
      if (0 != ffd_calculator(vector_function, x, ffd_jacobian_matrix))
        return -1;
      time_t ffd_end = clock();
      std::cout<<"FFD took "
               <<Scalar(ffd_end - ffd_begin) /
                 Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";
    }

    time_t ba_ffd_begin = clock();
    if (0 != ba_ffd_calculator(vector_function, x, ba_ffd_jacobian_matrix))
      return -1;
    time_t ba_ffd_end = clock();
    std::cout<<"BA FFD took "
             <<Scalar(ba_ffd_end - ba_ffd_begin) /
               Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";

    time_t ba_analytical_begin = clock();
    if (0 != ba_analytical_calculator(vector_function, x,
                                      ba_analytical_jacobian_matrix))
      return -1;
    time_t ba_analytical_end = clock();
    std::cout<<"BA Analytic took "
             <<Scalar(ba_analytical_end - ba_analytical_begin) /
               Scalar(CLOCKS_PER_SEC / 1000)<<" ms.\n";

    if (ba_analytical_jacobian_matrix.number_of_images() !=
        ba_ffd_jacobian_matrix.number_of_images() ||
        ba_analytical_jacobian_matrix.number_of_points() !=
        ba_ffd_jacobian_matrix.number_of_points() ||
        ba_analytical_jacobian_matrix.number_of_cameras() !=
        ba_ffd_jacobian_matrix.number_of_cameras() ||
        ba_analytical_jacobian_matrix.image_derivatives().size() !=
        ba_ffd_jacobian_matrix.image_derivatives().size() ||
        ba_analytical_jacobian_matrix.point_derivatives().size() !=
        ba_ffd_jacobian_matrix.point_derivatives().size() ||
        ba_analytical_jacobian_matrix.camera_derivatives().size() !=
        ba_ffd_jacobian_matrix.camera_derivatives().size())
    {
      std::cout<<"Number of images, points and cameras not match.\n";
      return -1;
    }

    Err result = 0;
    const Scalar threshold = Scalar(1e-3);
    for (Index i = 0;
         i < Index(ba_ffd_jacobian_matrix.point_derivatives().size());
         i++)
    {
      const typename BAAnalitycalJacobianMatrix::PointDerivativeBlock&
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
        std::cout<<"Size of point derivativematrix not match.\n";
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
          Scalar error = std::abs(ba_ffd_value - ba_analytical_value);
          if (ba_ffd_value != Scalar(0)) error = std::abs(error / ba_ffd_value);
          if (error > threshold)
          {
            std::cout<<"Deference beteen the point derivative matrix value"
                       "is too big\n";
            std::cout<<"error:"<<error<<"\n";
            std::cout<<"point derivative matrix id:"<<i<<".\n";
            std::cout<<"       ba_ffd_point_block["<<j<<","<<k<<"] = "
                     <<ba_ffd_value<<".\n";
            std::cout<<"ba_analytical_point_block["<<j<<","<<k<<"] = "
                     <<ba_analytical_value<<".\n";
            result = -1;
          }
        }
      }// for (Index j = 0; j < ba_ffd_point_block.derivative_block.rows(); j++)
    }

    for (Index i = 0;
         i < Index(ba_ffd_jacobian_matrix.image_derivatives().size());
         i++)
    {
      const typename BAAnalitycalJacobianMatrix::ImageDerivativeBlock&
        ba_analytical_image_block =
          ba_analytical_jacobian_matrix.image_derivatives()[i];
      const typename BAFFDJacobianMatrix::ImageDerivativeBlock&
        ba_ffd_image_block =
          ba_ffd_jacobian_matrix.image_derivatives()[i];

      if (ba_analytical_image_block.derivative_block.rows() !=
          ba_ffd_image_block.derivative_block.rows() ||
          ba_analytical_image_block.derivative_block.cols() !=
          ba_ffd_image_block.derivative_block.cols())
      {
        std::cout<<"Size of image derivativematrix not match.\n";
        result = -1;
        break;
      }

      for (Index j = 0; j < ba_ffd_image_block.derivative_block.rows(); j++)
      {
        for (Index k = 0; k < ba_ffd_image_block.derivative_block.cols(); k++)
        {
          Scalar ba_ffd_value =
            ba_ffd_image_block.derivative_block(j, k);
          Scalar ba_analytical_value =
            ba_analytical_image_block.derivative_block(j, k);
          Scalar error = std::abs(ba_ffd_value - ba_analytical_value);
          if (ba_ffd_value != Scalar(0)) error = std::abs(error / ba_ffd_value);
          if (error > threshold)
          {
            std::cout<<"Deference beteen the image derivative matrix value"
                       "is too big\n";
            std::cout<<"error:"<<error<<"\n";
            std::cout<<"image derivative matrix id:"<<i<<".\n";
            std::cout<<"       ba_ffd_image_block["<<j<<","<<k<<"] = "
                     <<ba_ffd_value<<".\n";
            std::cout<<"ba_analytical_image_block["<<j<<","<<k<<"] = "
                     <<ba_analytical_value<<".\n";
            result = -1;
          }
        }
      }// for (Index j = 0; j < ba_ffd_image_block.derivative_block.rows(); j++)
    }

    //for (Index i = 0;
    //     i < Index(ba_ffd_jacobian_matrix.camera_derivatives().size());
    //     i++)
    //{
    //  const typename BAAnalitycalJacobianMatrix::CameraDerivativeBlock&
    //    ba_analytical_camera_block =
    //      ba_analytical_jacobian_matrix.camera_derivatives()[i];
    //  const typename BAFFDJacobianMatrix::CameraDerivativeBlock&
    //    ba_ffd_camera_block =
    //      ba_ffd_jacobian_matrix.camera_derivatives()[i];

    //  if (ba_analytical_camera_block.derivative_block.rows() !=
    //      ba_ffd_camera_block.derivative_block.rows() ||
    //      ba_analytical_camera_block.derivative_block.cols() !=
    //      ba_ffd_camera_block.derivative_block.cols())
    //  {
    //    std::cout<<"Size of camera derivativematrix not match.\n";
    //    result = -1;
    //    break;
    //  }

    //  for (Index j = 0; j < ba_ffd_camera_block.derivative_block.rows(); j++)
    //  {
    //    for (Index k = 0; k < ba_ffd_camera_block.derivative_block.cols(); k++)
    //    {
    //      Scalar ba_ffd_value =
    //        ba_ffd_camera_block.derivative_block(j, k);
    //      Scalar ba_analytical_value =
    //        ba_analytical_camera_block.derivative_block(j, k);
    //      Scalar error = std::abs(ba_ffd_value - ba_analytical_value);
    //      if (ba_ffd_value != Scalar(0)) error = std::abs(error / ba_ffd_value);
    //      if (error > threshold)
    //      {
    //        std::cout<<"Deference beteen the camera derivative matrix value"
    //                   "is too big\n";
    //        std::cout<<"error:"<<error<<"\n";
    //        std::cout<<"camera derivative matrix id:"<<i<<".\n";
    //        std::cout<<"       ba_ffd_camera_block["<<j<<","<<k<<"] = "
    //                 <<ba_ffd_value<<".\n";
    //        std::cout<<"ba_analytical_camera_block["<<j<<","<<k<<"] = "
    //                 <<ba_analytical_value<<".\n";
    //        result = -1;
    //      }
    //    }
    //  }//for (Index j = 0; j < ba_ffd_camera_block.derivative_block.rows(); j++)
    //}

    if (with_ffd)
    {
      for (Index i = 0; i < y_size; i++)
      {
        for (Index j = 0; j < x_size; j++)
        {
          Scalar ffd_value = ffd_jacobian_matrix.coeff(i, j);
          Scalar ba_ffd_value = ba_ffd_jacobian_matrix.coeff(i, j);
          Scalar abs_error = std::abs(ba_ffd_value - ffd_value);
          if (abs_error > Scalar(1e-8))
          {
            std::cout<<"ba_ffd_jacobian_matrix["<<i<<","<<j<<"]="
                     <<ba_ffd_value<<".\n";
            std::cout<<"   ffd_jacobian_matrix["<<i<<","<<j<<"]="
                     <<ffd_value<<".\n";
            return -1;
          }
        }
      }
    }


    return 0;
  }
};

TEST(TestCameraSharedJacobianMatrixCalculator, SimpleTest)
{
  typedef double Scalar;
  typedef size_t ImageDimension;

  typedef TestCameraSharedJacobianMatrixCalculator<Scalar, ImageDimension>
          Tester;

  typedef hs::sfm::ba::CameraSharedSyntheticDataGenerator<Scalar,
                                                          ImageDimension>
          DataGenerator;
  typedef DataGenerator::FlightGenerator FlightGenerator;
  typedef DataGenerator::FlightGeneratorContainer FlightGeneratorContainer;
  typedef DataGenerator::IntrinsicParams IntrinsicParams;
  typedef DataGenerator::IntrinsicParamsContainer IntrinsicParamsContainer;

  typedef hs::sfm::ba::CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef VectorFunction::Index Index;
  typedef VectorFunction::XVector XVector;
  typedef VectorFunction::YVector YVector;

  FlightGeneratorContainer flight_generators;
  IntrinsicParamsContainer intrinsic_params_set;

  Scalar focal_length_in_metre_0 = 0.018858358970276164;
  size_t number_of_strips_0 = 3;
  size_t number_of_cameras_in_strip_0 = 3;
  Scalar ground_resolution_0 = 0.1;
  ImageDimension image_width_0 = 6000;
  ImageDimension image_height_0 = 4000;
  Scalar pixel_size_0 = 0.0000039;
  size_t number_of_points_0 = 2000;
  Scalar lateral_overlap_ratio_0 = 0.6;
  Scalar longitudinal_overlap_ratio_0 = 0.8;
  Scalar scene_max_height_0 = 50;
  Scalar camera_height_stddev_0 = 2;
  Scalar camera_planar_stddev_0 = 2;
  Scalar camera_rotation_stddev_0 = 10;
  FlightGenerator flight_generator_0(
    focal_length_in_metre_0,
    number_of_strips_0,
    number_of_cameras_in_strip_0,
    ground_resolution_0,
    image_width_0,
    image_height_0,
    pixel_size_0,
    number_of_points_0,
    lateral_overlap_ratio_0,
    longitudinal_overlap_ratio_0,
    scene_max_height_0,
    camera_height_stddev_0,
    camera_planar_stddev_0,
    camera_rotation_stddev_0);
  flight_generators.push_back(flight_generator_0);
  IntrinsicParams intrinsic_params_0(focal_length_in_metre_0 / pixel_size_0,
                                     0,
                                     -42.4095312016,
                                     -31.699212823,
                                     1,
                                     -0.0050490462006048831,
                                     0.031293804298609881,
                                     -0.030794820960442223,
                                     -0.00055376548320189127,
                                     -0.00049877717768381476);
  intrinsic_params_set.push_back(intrinsic_params_0);

  Scalar focal_length_in_metre_1 = 0.02995452167701055;
  size_t number_of_strips_1 = 4;
  size_t number_of_cameras_in_strip_1 = 3;
  Scalar ground_resolution_1 = 0.1;
  ImageDimension image_width_1 = 6000;
  ImageDimension image_height_1 = 4000;
  Scalar pixel_size_1 = 0.0000039;
  size_t number_of_points_1 = 2000;
  Scalar lateral_overlap_ratio_1 = 0.6;
  Scalar longitudinal_overlap_ratio_1 = 0.8;
  Scalar scene_max_height_1 = 50;
  Scalar camera_height_stddev_1 = 2;
  Scalar camera_planar_stddev_1 = 2;
  Scalar camera_rotation_stddev_1 = 10;
  FlightGenerator flight_generator_1(
    focal_length_in_metre_1,
    number_of_strips_1,
    number_of_cameras_in_strip_1,
    ground_resolution_1,
    image_width_1,
    image_height_1,
    pixel_size_1,
    number_of_points_1,
    lateral_overlap_ratio_1,
    longitudinal_overlap_ratio_1,
    scene_max_height_1,
    camera_height_stddev_1,
    camera_planar_stddev_1,
    camera_rotation_stddev_1);
  flight_generators.push_back(flight_generator_1);
  IntrinsicParams intrinsic_params_1(focal_length_in_metre_1 / pixel_size_1,
                                     0,
                                     -21.669436058,
                                     -44.8644764322,
                                     1,
                                     -0.02529179096221609,
                                     0.23762413973445157,
                                     -0.64208397668697237,
                                     -0.0020605099808780948,
                                     -0.00028706423764766859);
  intrinsic_params_set.push_back(intrinsic_params_1);

  Scalar focal_length_in_metre_2 = 0.019056097774998712;
  size_t number_of_strips_2 = 2;
  size_t number_of_cameras_in_strip_2 = 5;
  Scalar ground_resolution_2 = 0.1;
  ImageDimension image_width_2 = 6000;
  ImageDimension image_height_2 = 4000;
  Scalar pixel_size_2 = 0.0000039;
  size_t number_of_points_2 = 2000;
  Scalar lateral_overlap_ratio_2 = 0.6;
  Scalar longitudinal_overlap_ratio_2 = 0.8;
  Scalar scene_max_height_2 = 50;
  Scalar camera_height_stddev_2 = 2;
  Scalar camera_planar_stddev_2 = 2;
  Scalar camera_rotation_stddev_2 = 10;
  FlightGenerator flight_generator_2(
    focal_length_in_metre_2,
    number_of_strips_2,
    number_of_cameras_in_strip_2,
    ground_resolution_2,
    image_width_2,
    image_height_2,
    pixel_size_2,
    number_of_points_2,
    lateral_overlap_ratio_2,
    longitudinal_overlap_ratio_2,
    scene_max_height_2,
    camera_height_stddev_2,
    camera_planar_stddev_2,
    camera_rotation_stddev_2);
  flight_generators.push_back(flight_generator_2);
  IntrinsicParams intrinsic_params_2(focal_length_in_metre_2 / pixel_size_2,
                                     0,
                                     -35.2052431556,
                                     -16.4262220759,
                                     1,
                                     -0.10316088386868619,
                                     0.13520490482776426,
                                     -0.05489235547426094,
                                     4.1434720317373253e-006,
                                     -0.00025018439997095336);
  intrinsic_params_set.push_back(intrinsic_params_2);

  Scalar flight_longitudinal_overlap_ratio = 0.8;
  Scalar flight_lateral_overlap_ratio = 0.2;
  Scalar north_west_angle = 60.0;
  Scalar north_west_angle_stddev = 10.0;
  Scalar offset_stddev = 15.0;
  size_t number_of_points = 100;
  size_t number_of_planar_constrained_points = 4;
  size_t number_of_full_constrained_points = 6;
  size_t number_of_constrained_images = 10;
  size_t number_of_constrained_cameras = 2;

  DataGenerator data_generator(flight_longitudinal_overlap_ratio,
                               flight_lateral_overlap_ratio,
                               north_west_angle,
                               north_west_angle_stddev,
                               offset_stddev,
                               flight_generators,
                               number_of_points,
                               number_of_planar_constrained_points,
                               number_of_full_constrained_points,
                               number_of_constrained_images,
                               number_of_constrained_cameras,
                               intrinsic_params_set);

  VectorFunction vector_function;
  XVector x;
  YVector y_synthetic, y;

  ASSERT_EQ(0, data_generator(vector_function, x, y_synthetic));
  ASSERT_EQ(0, vector_function(x, y));

  Tester tester;
  ASSERT_EQ(0, tester.Test(vector_function, x, true));
}

}
