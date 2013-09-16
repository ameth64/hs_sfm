#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_FORWARD_FINITE_DIFFERENCE_JACOBIAN_MATRIX_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_FORWARD_FINITE_DIFFERENCE_JACOBIAN_MATRIX_CALCULATOR_HPP_

#include <cmath>

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_jacobian_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VectorFunction>
class BANaiveForwardFiniteDifferenceJacobianMatrixCalculator;

/**
 *  利用Forward Finit Differace计算简单的Bundle Adjustment函数的Jacobian矩阵
 */
template <typename _Scalar>
class BANaiveForwardFiniteDifferenceJacobianMatrixCalculator<
  BANaiveVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;

  typedef int Err;

  typedef BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::FeatureMap FeatureMap;
  typedef typename VectorFunction::FeatureMapContainer FeatureMapContainer;
  typedef BANaiveJacobianMatrix<Scalar, Index,
                                VectorFunction::params_per_feature_,
                                VectorFunction::params_per_camera_,
                                VectorFunction::params_per_point_>
          JacobianMatrix;
  typedef typename JacobianMatrix::DerivativeId DerivativeId;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;

  BANaiveForwardFiniteDifferenceJacobianMatrixCalculator()
    : delta_(Scalar(1e-4)), min_delta_(Scalar(1e-6)) {}

  BANaiveForwardFiniteDifferenceJacobianMatrixCalculator(
    Scalar delta, Scalar min_delta)
    : delta_(delta), min_delta_(min_delta) {}

  Err operator()(const VectorFunction& vector_function,
                 const XVector& x,
                 JacobianMatrix& jacobian_matrix) const
  {
    typedef Eigen::Triplet<DerivativeId, Index> TripletType;
    std::vector<TripletType> cameras_triplets;
    std::vector<TripletType> points_triplets;
    Index number_of_cameras = vector_function.number_of_cameras();
    Index number_of_points = vector_function.number_of_points();
    Index number_of_features = vector_function.number_of_features();
    Index camera_params_size = vector_function.GetCameraParamsSize();
    const FeatureMapContainer& feature_maps = vector_function.feature_maps();
    Index x_size = vector_function.GetXSize();
    jacobian_matrix.clear();
    jacobian_matrix.set_number_of_cameras(number_of_cameras);
    jacobian_matrix.set_number_of_points(number_of_points);
    for (Index i = 0; i < number_of_features; i++)
    {
      const FeatureMap& feature_map = feature_maps[i];
      Index j = feature_map.first;
      Index k = feature_map.second;

      Vector3 r = x.segment(j * VectorFunction::params_per_camera_, 3);
      Vector3 t = x.segment(j * VectorFunction::params_per_camera_ + 3, 3);
      Vector3 p = x.segment(camera_params_size + 
                            k * VectorFunction::params_per_point_,
                            VectorFunction::params_per_point_);
      Vector2 feature;
      VectorFunction::PointProjectToFeature(r, t, p, feature);

      typename JacobianMatrix::CameraDerivativeBlock camera_block;
      camera_block.camera_id = j;
      camera_block.feature_id = i;
      Vector3 r_delta = r;
      for (Index m = 0; m < 3; m++)
      {
        Scalar d = std::max(std::abs(delta_ * r[m]),
                            min_delta_);
        r_delta[m] += d;
        Vector2 feature_delta;
        VectorFunction::PointProjectToFeature(r_delta, t, p, feature_delta);
        for (Index n = 0; n < VectorFunction::params_per_feature_; n++)
        {
          camera_block.derivative_block.col(m)[n] =
            (feature_delta[n] - feature[n]) / d;
        }
        r_delta[m] = r[m];
      }
      Vector3 t_delta = t;
      for (Index m = 0; m < 3; m++)
      {
        Scalar d = std::max(std::abs(delta_ * t[m]),
                            min_delta_);
        t_delta[m] += d;
        Vector2 feature_delta;
        VectorFunction::PointProjectToFeature(r, t_delta, p, feature_delta);
        for (Index n = 0; n < VectorFunction::params_per_feature_; n++)
        {
          camera_block.derivative_block.col(3 + m)[n] =
            (feature_delta[n] - feature[n]) / d;
        }
        t_delta[m] = t[m];
      }
      jacobian_matrix.camera_derivatives().push_back(camera_block);
      cameras_triplets.push_back(
        TripletType(j, k, jacobian_matrix.camera_derivatives().size()));

      typename JacobianMatrix::PointDerivativeBlock point_block;
      point_block.point_id = k;
      point_block.feature_id = i;
      Vector3 p_delta = p;
      for (Index m = 0; m < VectorFunction::params_per_point_; m++)
      {
        Scalar d = std::max(std::abs(delta_ * p[m]),
                            min_delta_);
        p_delta[m] += d;
        Vector2 feature_delta;
        VectorFunction::PointProjectToFeature(r, t, p_delta, feature_delta);
        for (Index n = 0; n < VectorFunction::params_per_feature_; n++)
        {
          point_block.derivative_block.col(m)[n] =
            (feature_delta[n] - feature[n]) / d;
        }
        p_delta[m] = p[m];
      }
      jacobian_matrix.point_derivatives().push_back(point_block);
      points_triplets.push_back(
        TripletType(j, k, jacobian_matrix.point_derivatives().size()));
    }
    jacobian_matrix.camera_derivatives_map().resize(
      number_of_cameras, number_of_points);
    jacobian_matrix.camera_derivatives_map().setFromTriplets(
      cameras_triplets.begin(), cameras_triplets.end());
    jacobian_matrix.point_derivatives_map().resize(
      number_of_cameras, number_of_points);
    jacobian_matrix.point_derivatives_map().setFromTriplets(
      points_triplets.begin(), points_triplets.end());

    return 0;
  }

private:
  Scalar delta_;
  Scalar min_delta_;
};

}
}
}

#endif
