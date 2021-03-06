﻿#ifndef _HS_SFM_SYNTHETIC_MULTIPLE_FLIGHT_GENERATOR_HPP_
#define _HS_SFM_SYNTHETIC_MULTIPLE_FLIGHT_GENERATOR_HPP_

#include <cstdlib>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_sfm/synthetic/flight_generator.hpp"

namespace hs
{
namespace sfm
{
namespace synthetic
{

template <typename _Scalar, typename _ImageDimension>
class MultipleFlightGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

  typedef FlightGenerator<Scalar, ImageDimension> FlightGenerator;
  typedef EIGEN_STD_VECTOR(FlightGenerator) FlightGeneratorContainer; //16字节对齐的Eigen库的vector, 与STL兼容
  typedef typename FlightGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename ExtrinsicParams::Position Position;
  typedef typename ExtrinsicParams::Rotation Rotation;
  typedef typename FlightGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename FlightGenerator::Image Image;
  typedef typename FlightGenerator::ImageContainer ImageContainer;
  typedef typename FlightGenerator::Point3D Point3D;
  typedef typename FlightGenerator::Point3DContainer Point3DContainer;

private:
  typedef EIGEN_MATRIX(Scalar, 3, 3) RMatrix;
  typedef EIGEN_VECTOR(Scalar, 2) Point2D;
  typedef EIGEN_STD_VECTOR(Point2D) Point2DContainer;

public:
  MultipleFlightGenerator(Scalar flight_longitudinal_overlap_ratio, //纵向重叠率
                          Scalar flight_lateral_overlap_ratio,		//横向重叠率
                          Scalar north_west_angle,	//东偏北方位角
                          Scalar north_west_angle_stddev,	//方位角标准差
                          Scalar offset_stddev,	//航线平面平移标准差
                          const FlightGeneratorContainer& flight_generators)
    : flight_generators_(flight_generators)
  {
    size_t number_of_flights = flight_generators_.size();
    north_west_angles_.resize(number_of_flights);
    offsets_.resize(number_of_flights);
    Point2D offset = Point2D::Zero();
    for (size_t i = 0; i < number_of_flights; i++)
    {
      hs::math::random::NormalRandomVar<Scalar, 1>::Generate(
        north_west_angle,
        north_west_angle_stddev,
        north_west_angles_[i]);
      hs::math::random::NormalRandomVar<Scalar, 1>::Generate(
        offset[0],
        offset_stddev,
        offsets_[i][0]);
      hs::math::random::NormalRandomVar<Scalar, 1>::Generate(
        offset[1],
        offset_stddev,
        offsets_[i][1]);

      if (i < number_of_flights - 1)
      {
        const FlightGenerator& flight_this = flight_generators_[i];
        const FlightGenerator& flight_next = flight_generators_[i + 1];

        Scalar flight_this_x_dimension;
        Scalar flight_this_y_dimension;
        Scalar flight_this_z_dimension;
        flight_this.GetDimensions(flight_this_x_dimension,
                                  flight_this_y_dimension,
                                  flight_this_z_dimension);

        Scalar flight_next_x_dimension;
        Scalar flight_next_y_dimension;
        Scalar flight_next_z_dimension;
        flight_next.GetDimensions(flight_next_x_dimension,
                                  flight_next_y_dimension,
                                  flight_next_z_dimension);

        Scalar lateral_offset =
          (flight_this_x_dimension + flight_next_x_dimension) / Scalar(2) -
          flight_this_x_dimension * flight_lateral_overlap_ratio;
        Scalar longitudinal_offset =
          (flight_this_y_dimension + flight_next_y_dimension) / Scalar(2) -
          flight_this_y_dimension * flight_longitudinal_overlap_ratio;
        if (std::rand() % 2)
        {
          longitudinal_offset *= -1;
        }

        offset[0] += lateral_offset;
        offset[1] += longitudinal_offset;
      }//if (i < number_of_flights - 1)
    }//for (size_t i = 0; i < number_of_flights; i++)
  }

  Err GeneratePoints(size_t number_of_points,
                     Point3DContainer& points) const
  {
    points.resize(number_of_points);
    size_t number_of_flights = flight_generators_.size();
    Scalar pi = Scalar(3.141592653);
    for (size_t i = 0; i < number_of_points; i++)
    {
      size_t flight_id = i % number_of_flights;
      Scalar x_dimension, y_dimension, z_dimension;
      flight_generators_[flight_id].GetDimensions(x_dimension,
                                                  y_dimension,
                                                  z_dimension);
      Point3D max_point;
      max_point << x_dimension * 0.5,
                   y_dimension * 0.5,
                   z_dimension;
      Point3D min_point;
      min_point << -x_dimension * 0.5,
                   -y_dimension * 0.5,
                   0;
      hs::math::random::UniformRandomVar<Scalar, 3>::Generate(
        min_point, max_point, points[i]);
      points[i][0] += offsets_[flight_id][0];
      points[i][1] += offsets_[flight_id][1];
      RMatrix rotation;
      Scalar north_west_angle =
        north_west_angles_[flight_id] / Scalar(180) * pi;
      rotation << std::cos(north_west_angle), -std::sin(north_west_angle), 0,
                  std::sin(north_west_angle),  std::cos(north_west_angle), 0,
                  0, 0, 1;
      points[i] = rotation * points[i];
    }

    return 0;
  }

  Err GenerateExtrinsicParamsContainer(
    size_t flight_id,
    ExtrinsicParamsContainer& extrinsic_params_set,
    ImageContainer& images) const	//该方法计算相机参数并映射至世界坐标系(绝对坐标)
  {
    size_t number_of_flights = flight_generators_.size();
    Scalar pi = Scalar(3.141592653);

    const FlightGenerator& flight = flight_generators_[flight_id];
    flight.GenerateCameras(extrinsic_params_set, images);

    RMatrix rotation;	//计算相机在世界坐标系下的位移与方位角
    Scalar north_west_angle = north_west_angles_[flight_id] / Scalar(180) * pi;
    rotation << std::cos(north_west_angle), -std::sin(north_west_angle), 0,
                std::sin(north_west_angle),  std::cos(north_west_angle), 0,
                0, 0, 1;
    for (size_t i = 0; i < extrinsic_params_set.size(); i++)	//对每个图像的相机矩阵应用上述世界坐标系变换
    {
      ExtrinsicParams& extrinsic_params = extrinsic_params_set[i];
      extrinsic_params.position()[0] += offsets_[flight_id][0];
      extrinsic_params.position()[1] += offsets_[flight_id][1];
      extrinsic_params.position() = rotation * extrinsic_params.position();
      extrinsic_params.rotation() =
        RMatrix(extrinsic_params.rotation()) * rotation.transpose();
    }

    return 0;
  }

  size_t GetNumberOfFlights() const
  {
    return flight_generators_.size();
  }

private:
  /**
   *  保存多个驾次生成器的容器。
   */
  FlightGeneratorContainer flight_generators_;

  /**
   *  保存每个驾次对应的北偏东旋转角。
   */
  std::vector<Scalar> north_west_angles_;

  /**
   *  保存每个驾次对应的平面平移。
   */
   Point2DContainer offsets_;
};

}
}
}

#endif
