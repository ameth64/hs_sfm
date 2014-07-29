#ifndef _HS_SFM_CALIBRATE_PLANAR_CALIBRATION_GENERATOR_HPP_
#define _HS_SFM_CALIBRATE_PLANAR_CALIBRATION_GENERATOR_HPP_

#include <cmath>

#include "hs_math/random/uniform_random_var.hpp"
#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_math/geometry/euler_angles.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/projective_functions.hpp"
#include "hs_sfm/calibrate/planar_calibrator.hpp"

namespace hs
{
namespace sfm
{
namespace calibrate
{

/**
 *  平面标定模拟数据生成器。
 *
 *  该Functor生成模拟的平面点和影像点对应，以及多张相片拍摄时相机的外方位元素(
 *  即相机旋转与相机平移)，并为平面点与影像点加入高斯误差。其中平面点是位于z=0
 *  平面的棋盘点。相机外方位元素由以下方法生成：随机选取一个由三个旋转角组成的
 *  旋转，在该旋转后的坐标系的Z轴上取距离为f * w_p / w_i的点作为相机位置，其中
 *  f为焦距w_p为标定棋盘宽度，w_i为影像宽度，然后为相机位置以及相机旋转加入高斯
 *  误差以形成较为无规率的配置。
 */
template <typename _Scalar, typename _ImageDimension>
class PlanarCalibrationGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef hs::sfm::calibrate::PlanarCalibrator<Scalar> Calibrator;
  typedef typename Calibrator::IntrinsicParams IntrinsicParams;
  typedef typename Calibrator::ExtrinsicParams ExtrinsicParams;
  typedef typename Calibrator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
private:
  typedef typename Calibrator::Key Key;
  typedef typename Calibrator::Point Point;
  typedef typename Calibrator::Correspondence Correspondence;
  typedef typename Calibrator::PatternView PatternView;
  typedef typename Calibrator::PatternViewContainer PatternViewContainer;
  typedef hs::math::geometry::EulerAngles<Scalar> EulerAngles;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef typename EulerAngles::OrthoRotMat RMatrix;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

public:
  PlanarCalibrationGenerator(
    const IntrinsicParams& intrinsic_params,
    Scalar pattern_grid_size,
    size_t number_of_grid_rows,
    size_t number_of_grid_cols,
    size_t number_of_views,
    ImageDimension image_width,
    ImageDimension image_height)
    : intrinsic_params_(intrinsic_params),
      pattern_grid_size_(pattern_grid_size),
      number_of_grid_rows_(number_of_grid_rows),
      number_of_grid_cols_(number_of_grid_cols),
      number_of_views_(number_of_views),
      image_width_(image_width),
      image_height_(image_height) {}

  int operator() (ExtrinsicParamsContainer& extrinsic_params_set,
                  PatternViewContainer& pattern_views) const
  {
    extrinsic_params_set.resize(number_of_views_);
    pattern_views.resize(number_of_views_);

    Scalar pattern_width = pattern_grid_size_ * (number_of_grid_cols_ - 1);
    Scalar pattern_height = pattern_grid_size_ * (number_of_grid_rows_ - 1);
    Scalar distance =
      intrinsic_params_.focal_length() *
      std::max(pattern_width / Scalar(image_width_),
               pattern_height / Scalar(image_height_));
    Scalar camera_rotation_stddev = Scalar(5) * Scalar(M_PI) / Scalar(180);
    Scalar camera_position_x_stddev = distance * 0.02;
    Scalar camera_position_y_stddev = distance * 0.02;
    Scalar camera_position_z_stddev = distance * 0.04;
    for (size_t i = 0; i < number_of_views_; i++)
    {
      //生成相机旋转
      Vector3 euler_angles_mean, euler_angles_random, min, max;
      min.setZero();
      max << Scalar(20) * Scalar(M_PI) / Scalar(180),
             Scalar(20) * Scalar(M_PI) / Scalar(180),
             Scalar(M_PI);
      hs::math::random::UniformRandomVar<Scalar, 3>::Generate(
        min, max, euler_angles_mean);
      Matrix33 euler_angles_covariance = Matrix33::Identity();
      euler_angles_covariance *= camera_rotation_stddev *
                                 camera_rotation_stddev;
      hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
        euler_angles_mean, euler_angles_covariance, euler_angles_random);
      EulerAngles euler_angles;
      euler_angles[0] = euler_angles_random[0];
      euler_angles[1] = euler_angles_random[1];
      euler_angles[2] = euler_angles_random[2];
      RMatrix r_matrix = euler_angles.template ToOrthoRotMat<1, 2, 3, 1>();
      extrinsic_params_set[i].rotation() = r_matrix;

      //生成相机位置
      Vector3 camera_position_mean = r_matrix.col(2) * distance;
      Matrix33 camera_position_covariance = Matrix33::Identity();
      camera_position_covariance(0, 0) =
        camera_position_x_stddev * camera_position_x_stddev;
      camera_position_covariance(1, 1) =
        camera_position_y_stddev * camera_position_y_stddev;
      camera_position_covariance(2, 2) =
        camera_position_z_stddev * camera_position_z_stddev;
      hs::math::random::NormalRandomVar<Scalar, 3>::Generate(
        camera_position_mean, camera_position_covariance,
        extrinsic_params_set[i].position());

      std::cout<<"camera position:\n"<<extrinsic_params_set[i].position()<<"\n";

    }

    //生成平面三维点和影像二维点
    for (size_t i = 0; i < number_of_grid_rows_; i++)
    {
      for (size_t j = 0; j < number_of_grid_cols_; j++)
      {
        Point point;
        point << pattern_grid_size_ * Scalar(j) - pattern_width * Scalar(0.5),
                 pattern_grid_size_ * Scalar(i) - pattern_height * Scalar(0.5),
                 Scalar(0);

        for (size_t k = 0; k < number_of_views_; k++)
        {
          Key key =
            hs::sfm::ProjectiveFunctions<Scalar>::WorldPointProjectToImageKey(
              intrinsic_params_, extrinsic_params_set[k], point);
          if (key[0] >=0 && key[0] < Scalar(image_width_) &&
              key[1] >=0 && key[1] < Scalar(image_height_))
          {
            pattern_views[k].push_back(Correspondence(key, point));
          }
        }
      }
    }

    return 0;
  }

  const IntrinsicParams& intrinsic_params() const
  {
    return intrinsic_params_;
  }

  size_t number_of_views() const{return number_of_views_;}

private:
  /**
   *  相机内参数真实值
   */
  IntrinsicParams intrinsic_params_;
  /**
   *  标定棋盘一个格子的大小，以米为单位
   */
  Scalar pattern_grid_size_;
  /**
   *  标定棋盘格子行数
   */
  size_t number_of_grid_rows_;
  /**
   *  标定棋盘格子列数
   */
  size_t number_of_grid_cols_;
  /**
   *  随机生成的影像数量
   */
  size_t number_of_views_;

  /**
   *  影像宽度，像素为单位
   */
  ImageDimension image_width_;

  /**
   *  影像高度，像素为单位
   */
  ImageDimension image_height_;

};

}
}
}

#endif
