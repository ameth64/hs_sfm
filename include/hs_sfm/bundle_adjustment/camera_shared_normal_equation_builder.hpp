#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NORMAL_EQUATION_BUILDER_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NORMAL_EQUATION_BUILDER_HPP_

#include <vector>
#include <set>

#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_jacobian_matrix.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_normal_matrix.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_gradient.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_y_covariance_inverse.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedNormalEquationBuilder
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;
  typedef CameraSharedJacobianMatrix<Scalar> JacobianMatrix;
  typedef CameraSharedYCovarianceInverse<Scalar> YCovarianceInverse;
  typedef CameraSharedNormalMatrix<Scalar> NormalMatrix;
  typedef CameraSharedGradient<Scalar> Gradient;
  typedef YVector Residuals;

private:
  typedef typename JacobianMatrix::ImageDerivativeBlock ImageDerivativeBlock;
  typedef typename JacobianMatrix::PointDerivativeBlock PointDerivativeBlock;
  typedef typename JacobianMatrix::CameraDerivativeBlock CameraDerivativeBlock;
  typedef typename YCovarianceInverse::KeyBlock KeyBlock;
  typedef typename NormalMatrix::PointBlock PointBlock;
  typedef typename NormalMatrix::ImageBlock ImageBlock;
  typedef typename NormalMatrix::CameraBlock CameraBlock;
  typedef typename NormalMatrix::PointImageBlock PointImageBlock;
  typedef typename NormalMatrix::PointCameraBlock PointCameraBlock;
  typedef typename NormalMatrix::ImageCameraBlock ImageCameraBlock;
  typedef typename NormalMatrix::BlockIdx BlockIdx;
  class SparseTriplet
  {
  public:
    SparseTriplet() : row_(0), col_(0), value_(0) {}
    SparseTriplet(Index row, Index col, BlockIdx value)
      : row_(row), col_(col), value_(value) {}
    const Index& row() const {return row_;}
    const Index& col() const {return col_;}
    const BlockIdx& value() const {return value_;}
    bool operator == (const SparseTriplet& other) const
    {
      return row_ == other.row_ && col_ == other.col_;
    }
    bool operator != (const SparseTriplet& other) const
    {
      return !((*this) == other);
    }
    bool operator < (const SparseTriplet& other) const
    {
      return row_ < other.row_ || (!(row_ < other.row_) && col_ < other.col_);
    }
    bool operator <= (const SparseTriplet& other) const
    {
      return !(other < (*this));
    }
    bool operator > (const SparseTriplet& other) const
    {
      return (other < (*this));
    }
    bool operator >= (const SparseTriplet& other) const
    {
      return !((*this) < other);
    }

  private:
    Index row_;
    Index col_;
    BlockIdx value_;
  };
  typedef std::vector<SparseTriplet> SparseTripletVector;
  typedef std::set<SparseTriplet> SparseTripletSet;

public:
  Err operator() (const JacobianMatrix& jacobian_matrix,
                  const Residuals& residuals,
                  const YCovarianceInverse& y_covariance_inverse,
                  NormalMatrix& normal_matrix,
                  Gradient& gradient) const
  {
    if (ComputeNormalMatrix(jacobian_matrix,
                            y_covariance_inverse,
                            normal_matrix) != 0) return -1;

    if (ComputeGradient(jacobian_matrix,
                        residuals,
                        y_covariance_inverse,
                        gradient) != 0) return -1;

    return 0;
  }

public:
  Err ComputeNormalMatrix(const JacobianMatrix& jacobian_matrix,
                          const YCovarianceInverse& y_covariance_inverse,
                          NormalMatrix& normal_matrix) const
  {
    Index number_of_keys = jacobian_matrix.GetNumberOfKeys();
    size_t number_of_points = size_t(jacobian_matrix.number_of_points());
    size_t number_of_images = size_t(jacobian_matrix.number_of_images());
    size_t number_of_cameras = size_t(jacobian_matrix.number_of_cameras());
    Index camera_params_size =
      jacobian_matrix.GetIntrinsicParamsSizePerCamera();
    normal_matrix.Resize(number_of_points,
                         number_of_images,
                         number_of_cameras,
                         camera_params_size);
    SparseTripletVector point_image_triplets;
    SparseTripletSet point_camera_triplets_set;
    SparseTripletSet image_camera_triplets_set;
    //first pass, get triplets
    for (Index i = 0; i < number_of_keys; i++)
    {
      const PointDerivativeBlock& point_derivative =
        jacobian_matrix.point_derivatives()[i];
      const ImageDerivativeBlock& image_derivative =
        jacobian_matrix.image_derivatives()[i];
      const CameraDerivativeBlock& camera_derivative =
        jacobian_matrix.camera_derivatives()[i];

      Index point_id = point_derivative.point_id;
      Index image_id = image_derivative.image_id;
      Index camera_id = camera_derivative.camera_id;

      point_image_triplets.push_back(
        SparseTriplet(point_id, image_id, point_image_triplets.size() + 1));
      point_camera_triplets_set.insert(
        SparseTriplet(point_id, camera_id,
                      point_camera_triplets_set.size() + 1));
      image_camera_triplets_set.insert(
        SparseTriplet(image_id, camera_id,
                      image_camera_triplets_set.size() + 1));
    }
    //TODO:FIX ME!!Is copying elements to a vector nessasery?
    SparseTripletVector point_camera_triplets_vector(
      point_camera_triplets_set.size());
    std::copy(point_camera_triplets_set.begin(),
              point_camera_triplets_set.end(),
              point_camera_triplets_vector.begin());
    SparseTripletVector image_camera_triplets_vector(
      image_camera_triplets_set.size());
    std::copy(image_camera_triplets_set.begin(),
              image_camera_triplets_set.end(),
              image_camera_triplets_vector.begin());
    normal_matrix.SetPointImageMap(point_image_triplets.begin(),
                                   point_image_triplets.end());
    normal_matrix.SetPointCameraMap(point_camera_triplets_vector.begin(),
                                    point_camera_triplets_vector.end());
    normal_matrix.SetImageCameraMap(image_camera_triplets_vector.begin(),
                                    image_camera_triplets_vector.end());
    for (Index i = 0; i < number_of_keys; i++)
    {
      const KeyBlock& key_block = y_covariance_inverse.GetKeyBlock(i);

      const PointDerivativeBlock& point_derivative =
        jacobian_matrix.point_derivatives()[i];
      const ImageDerivativeBlock& image_derivative =
        jacobian_matrix.image_derivatives()[i];
      const CameraDerivativeBlock& camera_derivative =
        jacobian_matrix.camera_derivatives()[i];

      Index point_id = point_derivative.point_id;
      Index image_id = image_derivative.image_id;
      Index camera_id = camera_derivative.camera_id;

      PointBlock& point_block = normal_matrix.GetPointBlock(point_id);
      ImageBlock& image_block = normal_matrix.GetImageBlock(image_id);
      CameraBlock& camera_block = normal_matrix.GetCameraBlock(camera_id);

      point_block += point_derivative.derivative_block.transpose() *
                     key_block *
                     point_derivative.derivative_block;
      image_block += image_derivative.derivative_block.transpose() *
                     key_block *
                     image_derivative.derivative_block;
      camera_block += camera_derivative.derivative_block.transpose() *
                      key_block *
                      camera_derivative.derivative_block;

      PointImageBlock point_image_block =
        point_derivative.derivative_block.transpose() *
        key_block *
        image_derivative.derivative_block;
      if (normal_matrix.AddPointImageBlock(point_id, image_id,
                                           point_image_block) != 0) return -1;
      PointCameraBlock point_camera_block =
        point_derivative.derivative_block.transpose() *
        key_block *
        camera_derivative.derivative_block;
      if (normal_matrix.AddPointCameraBlock(point_id, camera_id,
                                            point_camera_block) != 0) return -1;
      ImageCameraBlock image_camera_block =
        image_derivative.derivative_block.transpose() *
        key_block *
        camera_derivative.derivative_block;
      if (normal_matrix.AddImageCameraBlock(image_id, camera_id,
                                            image_camera_block) != 0) return -1;
    }

    Index point_constaints_size = jacobian_matrix.GetPointConstraintsSize();
    for (Index offset = 0; offset < point_constaints_size; offset++)
    {
      Index point_id = jacobian_matrix.GetPointConstraintPointID(offset);
      Index param_id = jacobian_matrix.GetPointConstraintParamID(offset);
      normal_matrix.GetPointBlock(point_id)(param_id, param_id) +=
        y_covariance_inverse.GetConstraint(offset);
    }

    Index image_constaints_size = jacobian_matrix.GetImageConstraintsSize();
    for (Index offset = 0; offset < image_constaints_size; offset++)
    {
      Index image_id = jacobian_matrix.GetImageConstraintImageID(offset);
      Index param_id = jacobian_matrix.GetImageConstraintParamID(offset);
      normal_matrix.GetImageBlock(image_id)(param_id, param_id) +=
        y_covariance_inverse.GetConstraint(point_constaints_size + offset);
    }

    Index camera_constraints_size = jacobian_matrix.GetCameraConstraintsSize();
    for (Index offset = 0; offset < camera_constraints_size; offset++)
    {
      Index camera_id = jacobian_matrix.GetCameraConstraintCameraID(offset);
      Index param_id = jacobian_matrix.GetCameraConstraintParamID(offset);
      normal_matrix.GetCameraBlock(camera_id)(param_id, param_id) +=
        y_covariance_inverse.GetConstraint(point_constaints_size +
                                           image_constaints_size +
                                           offset);
    }

    return 0;
  }

  Err ComputeGradient(const JacobianMatrix& jacobian_matrix,
                      const Residuals& residuals,
                      const YCovarianceInverse& y_covariance_inverse,
                      Gradient& gradient) const
  {
    Index number_of_keys = jacobian_matrix.GetNumberOfKeys();
    size_t number_of_points = size_t(jacobian_matrix.number_of_points());
    size_t number_of_images = size_t(jacobian_matrix.number_of_images());
    size_t number_of_cameras = size_t(jacobian_matrix.number_of_cameras());
    Index camera_params_size =
      jacobian_matrix.GetIntrinsicParamsSizePerCamera();
    gradient.Resize(number_of_points,
                    number_of_images,
                    number_of_cameras,
                    camera_params_size);
    for (Index i = 0; i < number_of_keys; i++)
    {
      const KeyBlock& key_block = y_covariance_inverse.GetKeyBlock(i);

      const PointDerivativeBlock& point_derivative =
        jacobian_matrix.point_derivatives()[i];
      const ImageDerivativeBlock& image_derivative =
        jacobian_matrix.image_derivatives()[i];
      const CameraDerivativeBlock& camera_derivative =
        jacobian_matrix.camera_derivatives()[i];

      Index point_id = point_derivative.point_id;
      Index image_id = image_derivative.image_id;
      Index camera_id = camera_derivative.camera_id;

      gradient.GetPointSegment(point_id) +=
        point_derivative.derivative_block.transpose() *
        key_block *
        residuals.segment(i * VectorFunction::params_per_key_,
                          VectorFunction::params_per_key_);
      gradient.GetImageSegment(image_id) +=
        image_derivative.derivative_block.transpose() *
        key_block *
        residuals.segment(i * VectorFunction::params_per_key_,
                          VectorFunction::params_per_key_);
      gradient.GetCameraSegment(camera_id) +=
        camera_derivative.derivative_block.transpose() *
        key_block *
        residuals.segment(i * VectorFunction::params_per_key_,
                          VectorFunction::params_per_key_);
    }

    Index key_params_size = y_covariance_inverse.GetKeyParamsSize();
    Index point_constaints_size = jacobian_matrix.GetPointConstraintsSize();
    for (Index offset = 0; offset < point_constaints_size; offset++)
    {
      Index point_id = jacobian_matrix.GetPointConstraintPointID(offset);
      Index param_id = jacobian_matrix.GetPointConstraintParamID(offset);
      Index y_id = key_params_size + offset;
      gradient.GetPointSegment(point_id)[param_id] +=
        y_covariance_inverse.GetConstraint(offset) * residuals[y_id];
    }

    Index image_constaints_size = jacobian_matrix.GetImageConstraintsSize();
    for (Index offset = 0; offset < image_constaints_size; offset++)
    {
      Index image_id = jacobian_matrix.GetImageConstraintImageID(offset);
      Index param_id = jacobian_matrix.GetImageConstraintParamID(offset);
      Index constraint_id = point_constaints_size + offset;
      Index y_id = key_params_size + constraint_id;
      gradient.GetImageSegment(image_id)[param_id] +=
        y_covariance_inverse.GetConstraint(constraint_id) * residuals[y_id];
    }

    Index camera_constraints_size = jacobian_matrix.GetCameraConstraintsSize();
    for (Index offset = 0; offset < camera_constraints_size; offset++)
    {
      Index camera_id = jacobian_matrix.GetCameraConstraintCameraID(offset);
      Index param_id = jacobian_matrix.GetCameraConstraintParamID(offset);
      Index constraint_id =
        point_constaints_size + image_constaints_size + offset;
      Index y_id = key_params_size + constraint_id;
      gradient.GetCameraSegment(camera_id)[param_id] +=
        y_covariance_inverse.GetConstraint(constraint_id) * residuals[y_id];
    }

    return 0;
  }
};

}
}
}

#endif
