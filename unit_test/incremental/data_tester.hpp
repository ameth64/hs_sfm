#ifndef _HS_SFM_UNIT_TEST_INCREMENTAL_DATA_TESTER_HPP_
#define _HS_SFM_UNIT_TEST_INCREMENTAL_DATA_TESTER_HPP_

#include <iostream>
#include <fstream>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_math/geometry/euler_angles.hpp"

#include "hs_sfm/sfm_utility/key_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/incremental/reprojective_error_calculator.hpp"

namespace hs
{
namespace sfm
{
namespace incremental
{

template <typename _Scalar>
class DataTester
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;
  typedef hs::sfm::ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;

public:
  Err TestReprojectiveError(
    const KeysetContainer& keysets,
    const IntrinsicParamsContainer& intrinsic_params_set,
    const hs::sfm::TrackContainer& tracks,
    const hs::sfm::ObjectIndexMap& image_extrinsic_map,
    const hs::sfm::ObjectIndexMap& track_point_map,
    const hs::sfm::ViewInfoIndexer& view_info_indexer,
    const ExtrinsicParamsContainer& extrinsic_params_set,
    const Point3DContainer& points,
    Scalar key_stddev) const
  {
    typedef ReprojectiveErrorCalculator<Scalar> ReprojectiveErrorCalculator;
    ReprojectiveErrorCalculator reprojective_error_calculator;
    Scalar error = reprojective_error_calculator(
                     keysets,
                     intrinsic_params_set,
                     tracks,
                     image_extrinsic_map,
                     track_point_map,
                     view_info_indexer,
                     extrinsic_params_set,
                     points);
    std::cout<<"reprojective error:"<<error<<"\n";
    if (error < key_stddev + 1)
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }

  Err TestExtrinsicAccuracy(
    const ExtrinsicParamsContainer& extrinsic_params_set_true,
    const ExtrinsicParamsContainer& extrinsic_params_set_estimate,
    const std::string& accuracy_path,
    Scalar threshold) const
  {
    typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;
    size_t number_of_extrinsics = extrinsic_params_set_true.size();
    if (number_of_extrinsics != extrinsic_params_set_estimate.size())
    {
      return -1;
    }
    std::ofstream accuracy_file(accuracy_path.c_str(), std::ios::out);
    if (!accuracy_file.is_open())
    {
      return -1;
    }

    typedef hs::math::geometry::EulerAngles<Scalar> EulerAngles;
    Scalar mean_position_planar_error = Scalar(0);
    Scalar mean_position_height_error = Scalar(0);
    Scalar mean_rotation_angle0_error = Scalar(0);
    Scalar mean_rotation_angle1_error = Scalar(0);
    Scalar mean_rotation_angle2_error = Scalar(0);
    for (size_t i = 0; i < number_of_extrinsics; i++)
    {
      const ExtrinsicParams& extrinsic_params_estimate =
        extrinsic_params_set_estimate[i];
      const ExtrinsicParams& extrinsic_params_true =
        extrinsic_params_set_true[i];

      RMatrix rotation_estimate =
        extrinsic_params_estimate.rotation();
      EulerAngles angles_estimate;
      angles_estimate.template FromOrthoRotMat<2, 1, -3, 1>(
        rotation_estimate);

      RMatrix rotation_true =
        extrinsic_params_true.rotation();
      EulerAngles angles_true;
      angles_true.template FromOrthoRotMat<2, 1, -3, 1>(
        rotation_true);

      Scalar rotation_angle0_error =
          std::abs(angles_true[0] - angles_estimate[0]);
      Scalar rotation_angle1_error =
        std::abs(angles_true[1] - angles_estimate[1]);
      Scalar rotation_angle2_error =
        std::abs(angles_true[2] - angles_estimate[2]);

      const Point3D& position_estimate =
        extrinsic_params_estimate.position();
      const Point3D& position_true =
        extrinsic_params_true.position();

      Point3D position_diff = position_estimate - position_true;
      Scalar position_planar_error = position_diff.segment(0, 2).norm();
      Scalar position_height_error = std::abs(position_diff[2]);

      mean_position_planar_error += position_planar_error;
      mean_position_height_error += position_height_error;

      mean_rotation_angle0_error += rotation_angle0_error;
      mean_rotation_angle1_error += rotation_angle1_error;
      mean_rotation_angle2_error += rotation_angle2_error;

      accuracy_file<<i<<" "
                 <<position_estimate[0]<<" "
                 <<position_estimate[1]<<" "
                 <<position_estimate[2]<<" "
                 <<position_planar_error<<" "
                 <<position_height_error<<" "
                 <<rotation_angle0_error<<" "
                 <<rotation_angle1_error<<" "
                 <<rotation_angle2_error<<"\n";

      number_of_extrinsics++;
    }

    mean_position_planar_error /= Scalar(number_of_extrinsics);
    mean_position_height_error /= Scalar(number_of_extrinsics);
    mean_rotation_angle0_error /= Scalar(number_of_extrinsics);
    mean_rotation_angle1_error /= Scalar(number_of_extrinsics);
    mean_rotation_angle2_error /= Scalar(number_of_extrinsics);
    accuracy_file<<"mean position planar error:"
                 <<mean_position_planar_error<<"\n";
    accuracy_file<<"mean position height error:"
                 <<mean_position_height_error<<"\n";
    accuracy_file<<"mean rotation angle0 error:"
                 <<mean_rotation_angle0_error<<"\n";
    accuracy_file<<"mean rotation angle1 error:"
                 <<mean_rotation_angle1_error<<"\n";
    accuracy_file<<"mean rotation angle2 error:"
                 <<mean_rotation_angle2_error<<"\n";

    if (mean_position_height_error < threshold)
    {
      return 0;
    }
    else
    {
      return -1;
    }
  }

  Err TestPointsAccuracy(
    const Point3DContainer& points_true,
    const Point3DContainer& points_estimate,
    const std::string& accuracy_path,
    Scalar threshold) const
  {
    size_t number_of_points = points_true.size();
    if (points_estimate.size() != number_of_points)
    {
      return -1;
    }
    std::ofstream accuracy_file(accuracy_path.c_str(), std::ios::out);
    if (!accuracy_file.is_open())
    {
      return -1;
    }

    Scalar mean_point_planar_error = Scalar(0);
    Scalar mean_point_height_error = Scalar(0);
    for (size_t i = 0; i < number_of_points; i++)
    {
      const Point3D point_estimate = points_estimate[i];
      Point3D diff = point_estimate - points_true[i];
      Scalar planar_error = diff.segment(0, 2).norm();
      Scalar height_error = std::abs(diff(2));
      mean_point_planar_error += planar_error;
      mean_point_height_error += height_error;

      accuracy_file<<i<<" "
                   <<point_estimate[0]<<" "
                   <<point_estimate[1]<<" "
                   <<point_estimate[2]<<" "
                   <<planar_error<<" "
                   <<height_error<<"\n";
    }

    mean_point_planar_error /= Scalar(number_of_points);
    mean_point_planar_error /= Scalar(number_of_points);

    accuracy_file<<"mean point planar error:"
                 <<mean_point_planar_error<<"\n";
    accuracy_file<<"mean point height error:"
                 <<mean_point_height_error<<"\n";

    return 0;
  }
};

}
}
}

#endif
