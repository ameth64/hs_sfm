#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_FORWARD_FINITE_DIFFERENCE_JACOBIAN_MATRIX_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_FORWARD_FINITE_DIFFERENCE_JACOBIAN_MATRIX_CALCULATOR_HPP_

#include "hs_sfm/bundle_adjustment/camera_shared_jacobian_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VectorFunction>
class CameraSharedForwardFiniteDefferenceJacobianMatrixCalculator;

template <typename _Scalar>
class CameraSharedForwardFiniteDefferenceJacobianMatrixCalculator<
  CameraSharedVectorFunction<_Scalar> >
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;

  typedef CameraSharedJacobianMatrix<Scalar> JacobianMatrix;

private:
  typedef typename VectorFunction::ImageCameraMap ImageCameraMap;
  typedef typename VectorFunction::FeatureMap FeatureMap;
  typedef typename VectorFunction::FeatureMapContainer FeatureMapContainer;
  typedef typename JacobianMatrix::DerivativeId DerivativeId;
  typedef typename JacobianMatrix::KeyMap KeyMap;
  typedef typename JacobianMatrix::KeyMapContainer KeyMapContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) VectorX;

public:
  CameraSharedForwardFiniteDefferenceJacobianMatrixCalculator()
    : delta_(Scalar(1e-4)), min_delta_(Scalar(1e-6)) {}

  CameraSharedForwardFiniteDefferenceJacobianMatrixCalculator(
    Scalar delta, Scalar min_delta)
    : delta_(delta), min_delta_(min_delta) {}

  Err operator() (const VectorFunction& vector_function,
                  const XVector& x,
                  JacobianMatrix& jacobian_matrix) const
  {
    Index number_of_images = vector_function.number_of_images();
    Index number_of_points = vector_function.number_of_points();
    Index number_of_cameras = vector_function.number_of_cameras();
    Index number_of_keys = vector_function.number_of_keys();
    const FeatureMapContainer& feature_maps = vector_function.feature_maps();
    const ImageCameraMap& image_camera_map =
      vector_function.image_camera_map();
    Index x_size = vector_function.GetXSize();
    jacobian_matrix.Clear();
    jacobian_matrix.set_number_of_images(number_of_images);
    jacobian_matrix.set_number_of_points(number_of_points);
    jacobian_matrix.set_number_of_cameras(number_of_cameras);
    jacobian_matrix.intrinsic_computations_mask() =
      vector_function.intrinsic_computations_mask();
    Index point_params_size = vector_function.GetPointParamsSize();
    Index extrinsic_params_size = vector_function.GetExtrinsicParamsSize();
    Index intrinsic_params_size_per_camera =
      vector_function.GetIntrinsicParamsSizePerCamera();
    for (Index i = 0; i < number_of_keys; i++)
    {
      const FeatureMap& feature_map = feature_maps[i];
      Index image_id = feature_map.first;
      Index point_id = feature_map.second;
      Index camera_id = image_camera_map[image_id];
      KeyMap key_map;
      key_map.point_id = point_id;
      key_map.image_id = image_id;
      key_map.camera_id = camera_id;
      jacobian_matrix.key_maps().push_back(key_map);

      Vector3 point = x.segment(point_id * VectorFunction::params_per_point_,
                                VectorFunction::params_per_point_);
      Vector3 rotation =
        x.segment(point_params_size +
                  image_id * VectorFunction::extrinsic_params_per_image_, 3);
      Vector3 translation =
        x.segment(point_params_size +
                  image_id * VectorFunction::extrinsic_params_per_image_ + 3,
                  3);

      VectorX intrinsic_params =
        x.segment(point_params_size + extrinsic_params_size +
                  camera_id * intrinsic_params_size_per_camera,
                  intrinsic_params_size_per_camera);
      Vector2 image_key = vector_function.WorldPointToImageKey(
                            point, rotation, translation, intrinsic_params);

      typename JacobianMatrix::PointDerivativeBlock point_block;
      point_block.point_id = point_id;
      point_block.key_id = i;
      Vector3 point_delta = point;
      for (Index j = 0; j < VectorFunction::params_per_point_; j++)
      {
        Scalar d = std::max(std::abs(delta_ * point[j]),
                            min_delta_);
        point_delta[j] += d;
        Vector2 image_key_delta =
          vector_function.WorldPointToImageKey(point_delta,
                                               rotation,
                                               translation,
                                               intrinsic_params);
        for (Index k = 0; k < VectorFunction::params_per_key_; k++)
        {
          point_block.derivative_block.col(j)[k] =
            (image_key_delta[k] - image_key[k]) / d;
        }
        point_delta[j] = point[j];
      }
      jacobian_matrix.point_derivatives().push_back(point_block);

      typename JacobianMatrix::ImageDerivativeBlock image_block;
      image_block.image_id = image_id;
      image_block.key_id = i;
      Vector3 rotation_delta = rotation;
      for (Index j = 0; j < 3; j++)
      {
        Scalar d = std::max(std::abs(delta_ * rotation[j]),
                            min_delta_);
        rotation_delta[j] += d;
        Vector2 image_key_delta =
          vector_function.WorldPointToImageKey(point,
                                               rotation_delta,
                                               translation,
                                               intrinsic_params);
        for(Index k = 0; k < VectorFunction::params_per_key_; k++)
        {
          image_block.derivative_block.col(j)[k] =
            (image_key_delta[k] - image_key[k]) / d;
        }
        rotation_delta[j] = rotation[j];
      }
      Vector3 translation_delta = translation;
      for (Index j = 0; j < 3; j++)
      {
        Scalar d = std::max(std::abs(delta_ * translation[j]),
                            min_delta_);
        translation_delta[j] += d;
        Vector2 image_key_delta =
          vector_function.WorldPointToImageKey(point,
                                               rotation,
                                               translation_delta,
                                               intrinsic_params);
        for (Index k = 0; k < VectorFunction::params_per_key_; k++)
        {
          image_block.derivative_block.col(3 + j)[k] =
            (image_key_delta[k] - image_key[k]) / d;
        }
        translation_delta[j] = translation[j];
      }
      jacobian_matrix.image_derivatives().push_back(image_block);

      typename JacobianMatrix::CameraDerivativeBlock camera_block;
      camera_block.derivative_block.resize(VectorFunction::params_per_key_,
                                           intrinsic_params_size_per_camera);
      camera_block.camera_id = camera_id;
      camera_block.key_id = i;
      VectorX intrinsic_params_delta = intrinsic_params;
      for (Index j = 0; j < intrinsic_params_size_per_camera; j++)
      {
        Scalar d = std::max(std::abs(delta_ * intrinsic_params[j]),
                            min_delta_);
        intrinsic_params_delta[j] += d;
        Vector2 image_key_delta =
          vector_function.WorldPointToImageKey(point,
                                               rotation,
                                               translation,
                                               intrinsic_params_delta);
        for (Index k = 0; k < VectorFunction::params_per_key_; k++)
        {
          camera_block.derivative_block.col(j)[k] =
            (image_key_delta[k] - image_key[k]) / d;
        }
        intrinsic_params_delta[j] = intrinsic_params[j];
      }
      jacobian_matrix.camera_derivatives().push_back(camera_block);
    }// for (Index i = 0; i < number_of_keys; i++)

    jacobian_matrix.SetPointConstraints(vector_function.point_constraints());
    jacobian_matrix.SetImageConstraints(vector_function.image_constraints());
    jacobian_matrix.SetCameraConstraints(vector_function.camera_constraints());

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
