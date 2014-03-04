#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_JACOBIAN_MATRIX_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_JACOBIAN_MATRIX_HPP_

#include <limits>
#include <vector>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedJacobianMatrix
{
public:
  typedef _Scalar Scalar;

  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;

  struct ImageDerivativeBlock
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef EIGEN_MATRIX(Scalar,
                         VectorFunction::params_per_key_,
                         VectorFunction::extrinsic_params_per_image_)
            DerivativeBlock;
    Index image_id;
    Index key_id;
    DerivativeBlock derivative_block;
  };

  struct PointDerivativeBlock
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef EIGEN_MATRIX(Scalar,
                         VectorFunction::params_per_key_,
                         VectorFunction::params_per_point_)
            DerivativeBlock;
    Index point_id;
    Index key_id;
    DerivativeBlock derivative_block;
  };

  struct CameraDerivativeBlock
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef EIGEN_MATRIX(Scalar,
                         VectorFunction::params_per_key_,
                         Eigen::Dynamic)
            DerivativeBlock;
    Index camera_id;
    Index key_id;
    DerivativeBlock derivative_block;
  };

  typedef EIGEN_STD_VECTOR(ImageDerivativeBlock)
          ImageDerivativeBlockContainer;
  typedef EIGEN_STD_VECTOR(PointDerivativeBlock)
          PointDerivativeBlockContainer;
  typedef EIGEN_STD_VECTOR(CameraDerivativeBlock)
          CameraDerivativeBlockContainer;
  enum
  {
    EigenDefaultMajor =
#if defined(EIGEN_DEFAULT_TO_ROW_MAJOR)
    Eigen::RowMajor
#else
    Eigen::ColMajor
#endif
  };
  typedef typename ImageDerivativeBlockContainer::size_type DerivativeId;
  typedef EIGEN_SPARSE_MATRIX(DerivativeId, EigenDefaultMajor, Index)
          DerivativeMap;

  struct KeyMap
  {
    Index point_id;
    Index image_id;
    Index camera_id;
  };
  typedef std::vector<KeyMap> KeyMapContainer;

public:
  inline void Clear()
  {
    ImageDerivativeBlockContainer().swap(image_derivatives_);
    PointDerivativeBlockContainer().swap(point_derivatives_);
    CameraDerivativeBlockContainer().swap(camera_derivatives_);
    KeyMapContainer().swap(key_maps_);
  }

  Scalar coeff(Index i, Index j) const
  {
    Index keys_params_size = GetKeysParamsSize();
    Index point_constraints_end = keys_params_size + GetPointConstraintsSize();
    Index image_constraints_end = point_constraints_end +
                                  GetImageConstraintsSize();
    Index camera_constraints_end = image_constraints_end +
                                   GetCameraConstraintsSize();
    if (keys_params_size == 0)
    {
      return std::numeric_limits<Scalar>::signaling_NaN();
    }

    if (i < keys_params_size)
    {
      Index key_id = i / VectorFunction::params_per_key_;
      Index key_param_id = i % VectorFunction::params_per_key_;

      Index points_params_size = GetPointParamsSize();
      Index extrinsic_params_size = GetExtrinsicParamsSize();
      Index intrinsic_params_size = GetIntrinsicParamsSize();
      Index points_params_begin = 0;
      Index points_params_end =
        point_derivatives_.empty() ? 0 : points_params_size;
      Index images_params_begin = points_params_end;
      Index images_params_end =
        images_params_begin +
        (image_derivatives_.empty() ? 0 : extrinsic_params_size);
      Index cameras_params_begin = images_params_end;
      Index cameras_params_end =
        cameras_params_begin +
        (camera_derivatives_.empty() ? 0 : intrinsic_params_size);
      if ((!point_derivatives_.empty()) && j < points_params_end)
      {
        Index offset = j - points_params_begin;
        Index point_id = offset / VectorFunction::params_per_point_;
        Index point_param_id = offset % VectorFunction::params_per_point_;
        const PointDerivativeBlock& point_derivative_block =
          point_derivatives_[key_id];
        if (point_derivative_block.point_id == point_id)
        {
          return point_derivative_block.derivative_block(
                   key_param_id, point_param_id);
        }
        else
        {
          return Scalar(0);
        }
      }
      else if ((!image_derivatives_.empty()) && j < images_params_end)
      {
        Index offset = j - images_params_begin;
        Index image_id = offset / VectorFunction::extrinsic_params_per_image_;
        Index image_param_id =
          offset % VectorFunction::extrinsic_params_per_image_;
        const ImageDerivativeBlock& image_derivative_block =
          image_derivatives_[key_id];
        if (image_derivative_block.image_id == image_id)
        {
          return image_derivative_block.derivative_block(
                   key_param_id, image_param_id);
        }
        else
        {
          return Scalar(0);
        }
      }
      else if ((!camera_derivatives_.empty()) && j < cameras_params_end)
      {
        Index offset = j - cameras_params_begin;
        Index camera_id = offset / GetIntrinsicParamsSizePerCamera();
        Index camera_param_id = offset % GetIntrinsicParamsSizePerCamera();
        const CameraDerivativeBlock& camera_derivative_block =
          camera_derivatives_[key_id];
        if (camera_derivative_block.camera_id == camera_id)
        {
          return camera_derivative_block.derivative_block(
                   key_param_id, camera_param_id);
        }
        else
        {
          return Scalar(0);
        }
      }
      else
      {
        return std::numeric_limits<Scalar>::signaling_NaN();
      }
    }// if (i < keys_params_size)
    else if (i < point_constraints_end)
    {
      Index y_offset = i - keys_params_size;
      if (point_constraints_map_[y_offset] == j)
      {
        return Scalar(1);
      }
      else
      {
        return Scalar(0);
      }
    }// else if (i < point_constraints_end)
    else if (i < image_constraints_end)
    {
      Index y_offset = i - point_constraints_end;
      if (image_constraints_map_[y_offset] == j)
      {
        return Scalar(1);
      }
      else
      {
        return Scalar(0);
      }
    }// else if (i < image_constraints_end)
    else
    {
      Index y_offset = i - image_constraints_end;
      if (camera_constraints_map_[y_offset] == j)
      {
        return Scalar(1);
      }
      else
      {
        return Scalar(0);
      }
    }
  }

  inline Index GetExtrinsicParamsSize() const
  {
    return number_of_images_ * VectorFunction::extrinsic_params_per_image_;
  }

  inline Index GetPointParamsSize() const
  {
    return number_of_points_ * VectorFunction::params_per_point_;
  }

  inline Index GetIntrinsicParamsSizePerCamera() const
  {
    Index params_per_camera = 0;
    if (intrinsic_computations_mask_[COMPUTE_RADIAL_DISTORTION])
    {
      params_per_camera += 3;
    }
    if (intrinsic_computations_mask_[COMPUTE_DECENTERING_DISTORTION])
    {
      params_per_camera += 2;
    }
    if (intrinsic_computations_mask_[COMPUTE_INTRINSIC_PARAMS])
    {
      params_per_camera += 5;
    }
    return params_per_camera;
  }

  inline Index GetIntrinsicParamsSize() const
  {
    return number_of_cameras_ * GetIntrinsicParamsSizePerCamera();
  }

  inline Index GetNumberOfKeys() const
  {
    return Index(key_maps_.size());
  }

  inline Index GetKeysParamsSize() const
  {
    return GetNumberOfKeys() * VectorFunction::params_per_key_;
  }

  inline Index GetPointConstraintsSize() const
  {
    if (point_derivatives_.empty())
    {
      return 0;
    }
    else
    {
      return Index(point_constraints_map_.size());
    }
  }

  inline Index GetImageConstraintsSize() const
  {
    if (image_derivatives_.empty())
    {
      return 0;
    }
    else
    {
      return Index(image_constraints_map_.size());
    }
  }

  inline Index GetCameraConstraintsSize() const
  {
    if (camera_derivatives_.empty())
    {
      return 0;
    }
    else
    {
      return Index(camera_constraints_map_.size());
    }
  }

  inline Index GetXSize() const
  {
    return ((point_derivatives_.empty() ? 0 : GetPointParamsSize()) +
            (image_derivatives_.empty() ? 0 : GetExtrinsicParamsSize()) +
            (camera_derivatives_.empty() ? 0 : GetIntrinsicParamsSize()));
  }

  inline Index GetYSize() const
  {
    return GetKeysParamsSize() +
           GetPointConstraintsSize() +
           GetImageConstraintsSize() +
           GetCameraConstraintsSize();
  }

  Index number_of_images() const
  {
    return number_of_images_;
  }
  void set_number_of_images(Index number_of_images)
  {
    number_of_images_ = number_of_images;
  }

  Index number_of_points() const
  {
    return number_of_points_;
  }
  void set_number_of_points(Index number_of_points)
  {
    number_of_points_ = number_of_points;
  }

  Index number_of_cameras() const
  {
    return number_of_cameras_;
  }
  void set_number_of_cameras(Index number_of_cameras)
  {
    number_of_cameras_ = number_of_cameras;
  }

  const IntrinsicComputationsMask& intrinsic_computations_mask() const
  {
    return intrinsic_computations_mask_;
  }
  IntrinsicComputationsMask& intrinsic_computations_mask()
  {
    return intrinsic_computations_mask_;
  }

  const ImageDerivativeBlockContainer& image_derivatives() const
  {
    return image_derivatives_;
  }
  ImageDerivativeBlockContainer& image_derivatives()
  {
    return image_derivatives_;
  }

  const PointDerivativeBlockContainer& point_derivatives() const
  {
    return point_derivatives_;
  }
  PointDerivativeBlockContainer& point_derivatives()
  {
    return point_derivatives_;
  }

  const CameraDerivativeBlockContainer& camera_derivatives() const
  {
    return camera_derivatives_;
  }
  CameraDerivativeBlockContainer& camera_derivatives()
  {
    return camera_derivatives_;
  }

  const KeyMapContainer& key_maps() const
  {
    return key_maps_;
  }
  KeyMapContainer& key_maps()
  {
    return key_maps_;
  }

  void SetPointConstraints(const PointConstraintContainer& point_constraints)
  {
    point_constraints_map_.clear();
    auto itr_point = point_constraints.begin();
    auto itr_point_end = point_constraints.end();
    for (; itr_point != itr_point_end; ++itr_point)
    {
      Index point_id = Index(itr_point->point_id);
      Index x_offset = point_id * VectorFunction::params_per_point_;
      if (itr_point->mask[POINT_CONSTRAIN_X])
      {
        point_constraints_map_.push_back(x_offset + 0);
      }
      if (itr_point->mask[POINT_CONSTRAIN_Y])
      {
        point_constraints_map_.push_back(x_offset + 1);
      }
      if (itr_point->mask[POINT_CONSTRAIN_Z])
      {
        point_constraints_map_.push_back(x_offset + 2);
      }
    }
  }

  void SetImageConstraints(const ImageConstraintContainer& image_constraints)
  {
    image_constraints_map_.clear();
    auto itr_image = image_constraints.begin();
    auto itr_image_end = image_constraints.end();
    Index image_params_begin =
      point_derivatives_.empty() ? 0 : GetPointParamsSize();
    for (; itr_image != itr_image_end; ++itr_image)
    {
      Index image_id = Index(itr_image->image_id);
      Index x_offset = image_params_begin +
                       image_id * VectorFunction::extrinsic_params_per_image_;
      if (itr_image->mask[IMAGE_CONSTRAIN_ROTATION])
      {
        image_constraints_map_.push_back(x_offset + 0);
        image_constraints_map_.push_back(x_offset + 1);
        image_constraints_map_.push_back(x_offset + 2);
      }
      if (itr_image->mask[IMAGE_CONSTRAIN_POSITION_X])
      {
        image_constraints_map_.push_back(x_offset + 3);
      }
      if (itr_image->mask[IMAGE_CONSTRAIN_POSITION_Y])
      {
        image_constraints_map_.push_back(x_offset + 4);
      }
      if (itr_image->mask[IMAGE_CONSTRAIN_POSITION_Z])
      {
        image_constraints_map_.push_back(x_offset + 5);
      }
    }
  }

  void SetCameraConstraints(const CameraConstraintContainer& camera_constraints)
  {
    camera_constraints_map_.clear();
    auto itr_camera = camera_constraints.begin();
    auto itr_camera_end = camera_constraints.end();
    Index camera_params_begin =
      (point_derivatives_.empty() ? 0 : GetPointParamsSize()) +
      (image_derivatives_.empty() ? 0 : GetExtrinsicParamsSize());
    Index intrinsic_params_per_camera = GetIntrinsicParamsSizePerCamera();
    for (; itr_camera != itr_camera_end; ++itr_camera)
    {
      Index camera_id = Index(itr_camera->camera_id);
      Index x_offset = camera_params_begin +
                       camera_id * intrinsic_params_per_camera;
      if (itr_camera->radial_mask[RADIAL_CONSTRAIN_K1])
      {
        camera_constraints_map_.push_back(x_offset + 0);
      }
      if (itr_camera->radial_mask[RADIAL_CONSTRAIN_K2])
      {
        camera_constraints_map_.push_back(x_offset + 1);
      }
      if (itr_camera->radial_mask[RADIAL_CONSTRAIN_K3])
      {
        camera_constraints_map_.push_back(x_offset + 2);
      }
      if (intrinsic_computations_mask_[COMPUTE_RADIAL_DISTORTION])
      {
        x_offset += 3;
      }
      if (itr_camera->decentering_mask[DECENTERING_CONSTRAIN_D1])
      {
        camera_constraints_map_.push_back(x_offset + 0);
      }
      if (itr_camera->decentering_mask[DECENTERING_CONSTRAIN_D2])
      {
        camera_constraints_map_.push_back(x_offset + 1);
      }
      if (intrinsic_computations_mask_[COMPUTE_DECENTERING_DISTORTION])
      {
        x_offset += 2;
      }
      if (itr_camera->intrinsic_mask[INTRINSIC_CONSTRAIN_FOCAL_LENGTH])
      {
        camera_constraints_map_.push_back(x_offset + 0);
      }
      if (itr_camera->intrinsic_mask[INTRINSIC_CONSTRAIN_SKEW])
      {
        camera_constraints_map_.push_back(x_offset + 1);
      }
      if (itr_camera->intrinsic_mask[INTRINSIC_CONSTRAIN_PRINCIPAL_X])
      {
        camera_constraints_map_.push_back(x_offset + 2);
      }
      if (itr_camera->intrinsic_mask[INTRINSIC_CONSTRAIN_PRINCIPAL_Y])
      {
        camera_constraints_map_.push_back(x_offset + 3);
      }
      if (itr_camera->intrinsic_mask[INTRINSIC_CONSTRAIN_PIXEL_RATIO])
      {
        camera_constraints_map_.push_back(x_offset + 4);
      }
    }// for (; itr_camera != itr_camera_end; ++itr_camera)
  }

  Index GetPointConstraintPointID(Index offset) const
  {
    return point_constraints_map_[offset] /
           VectorFunction::params_per_point_;
  }
  Index GetPointConstraintParamID(Index offset) const
  {
    return point_constraints_map_[offset] %
           VectorFunction::params_per_point_;
  }

  Index GetImageConstraintImageID(Index offset) const
  {
    Index image_params_begin =
      point_derivatives_.empty() ? 0 : GetPointParamsSize();
    return (image_constraints_map_[offset] - image_params_begin) /
           VectorFunction::extrinsic_params_per_image_;
  }
  Index GetImageConstraintParamID(Index offset) const
  {
    Index image_params_begin =
      point_derivatives_.empty() ? 0 : GetPointParamsSize();
    return (image_constraints_map_[offset] - image_params_begin) %
           VectorFunction::extrinsic_params_per_image_;
  }

  Index GetCameraConstraintCameraID(Index offset) const
  {
    Index camera_params_begin =
      (point_derivatives_.empty() ? 0 : GetPointParamsSize()) +
      (image_derivatives_.empty() ? 0 : GetExtrinsicParamsSize());
    return (camera_constraints_map_[offset] - camera_params_begin) /
           GetIntrinsicParamsSizePerCamera();
  }
  Index GetCameraConstraintParamID(Index offset) const
  {
    Index camera_params_begin =
      (point_derivatives_.empty() ? 0 : GetPointParamsSize()) +
      (image_derivatives_.empty() ? 0 : GetExtrinsicParamsSize());
    return (camera_constraints_map_[offset] - camera_params_begin) %
           GetIntrinsicParamsSizePerCamera();
  }

private:
  ImageDerivativeBlockContainer image_derivatives_;
  PointDerivativeBlockContainer point_derivatives_;
  CameraDerivativeBlockContainer camera_derivatives_;
  KeyMapContainer key_maps_;
  Index number_of_images_;
  Index number_of_points_;
  Index number_of_cameras_;
  IntrinsicComputationsMask intrinsic_computations_mask_;
  std::vector<Index> point_constraints_map_;
  std::vector<Index> image_constraints_map_;
  std::vector<Index> camera_constraints_map_;
};

}
}
}

#endif
