#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NORMAL_MATRIX_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NORMAL_MATRIX_HPP_

#include <limits>
#include <algorithm>

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedNormalMatrix
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;

  static const Index params_per_key_ = VectorFunction::params_per_key_;
  static const Index params_per_point_ = VectorFunction::params_per_point_;
  static const Index extrinsic_params_per_image_ =
    VectorFunction::extrinsic_params_per_image_;

  typedef EIGEN_MATRIX(Scalar,
                       params_per_point_,
                       params_per_point_) PointBlock;
  typedef EIGEN_MATRIX(Scalar,
                       extrinsic_params_per_image_,
                       extrinsic_params_per_image_) ImageBlock;
  typedef EIGEN_MATRIX(Scalar,
                       Eigen::Dynamic,
                       Eigen::Dynamic) CameraBlock;
  typedef EIGEN_MATRIX(Scalar,
                       params_per_point_,
                       extrinsic_params_per_image_) PointImageBlock;
  typedef EIGEN_MATRIX(Scalar,
                       params_per_point_,
                       Eigen::Dynamic) PointCameraBlock;
  typedef EIGEN_MATRIX(Scalar,
                       extrinsic_params_per_image_,
                       Eigen::Dynamic) ImageCameraBlock;
  typedef EIGEN_STD_VECTOR(PointBlock) PointBlockContainer;
  typedef EIGEN_STD_VECTOR(ImageBlock) ImageBlockContainer;
  typedef EIGEN_STD_VECTOR(CameraBlock) CameraBlockContainer;
  typedef EIGEN_STD_VECTOR(PointImageBlock) PointImageBlockContainer;
  typedef EIGEN_STD_VECTOR(PointCameraBlock) PointCameraBlockContainer;
  typedef EIGEN_STD_VECTOR(ImageCameraBlock) ImageCameraBlockContainer;

  enum
  {
    EigenDefaultMajor =
#if defined(EIGEN_DEFAULT_TO_ROW_MAJOR)
      Eigen::RowMajor
#else
      Eigen::ColMajor
#endif
  };
  typedef typename PointBlockContainer::size_type BlockIdx;
  typedef EIGEN_SPARSE_MATRIX(BlockIdx, EigenDefaultMajor, Index)
          SparseBlockMap;

public:
  Scalar coeff(Index i, Index j) const
  {
    Index points_end = Index(number_of_points_) * params_per_point_;
    Index images_end =
      points_end + Index(number_of_images_) * extrinsic_params_per_image_;
    Index cameras_end =
      images_end + Index(number_of_cameras_) * camera_params_size_;

    if (i < 0 || i >= cameras_end ||
        j < 0 || j >= cameras_end)
    {
      return std::numeric_limits<Scalar>::signaling_NaN();
    }
    else if (i < points_end && j < points_end &&
             i / params_per_point_ == j / params_per_point_)
    {
      Index point_id = i / params_per_point_;
      Index point_row_offset = i % params_per_point_;
      Index point_col_offset = j % params_per_point_;
      return point_blocks_[point_id](point_row_offset, point_col_offset);
    }
    else if (i >= points_end && j >= points_end &&
             i < images_end && j < images_end &&
             (i - points_end) / extrinsic_params_per_image_ ==
             (j - points_end) / extrinsic_params_per_image_)
    {
      Index image_id = (i - points_end) / extrinsic_params_per_image_;
      Index image_row_offset = (i - points_end) % extrinsic_params_per_image_;
      Index image_col_offset = (j - points_end) % extrinsic_params_per_image_;
      return image_blocks_[image_id](image_row_offset, image_col_offset);
    }
    else if (i >= images_end && j >= images_end &&
             i < cameras_end && j < cameras_end &&
             (i - images_end) / camera_params_size_ ==
             (j - images_end) / camera_params_size_)
    {
      Index camera_id = (i - images_end) / camera_params_size_;
      Index camera_row_offset = (i - images_end) % camera_params_size_;
      Index camera_col_offset = (j - images_end) % camera_params_size_;
      return camera_blocks_[camera_id](camera_row_offset, camera_col_offset);
    }
    else if (i < points_end && j >= points_end && j < images_end)
    {
      Index point_id = i / params_per_point_;
      Index image_id = (j - points_end) / extrinsic_params_per_image_;
      Index row_offset = i % params_per_point_;
      Index col_offset = (j - points_end) % extrinsic_params_per_image_;
      const PointImageBlock* block_ptr = GetPointImageBlock(point_id,
                                                            image_id);
      if (block_ptr == NULL)
      {
        return Scalar(0);
      }
      else
      {
        return (*block_ptr)(row_offset, col_offset);
      }
    }
    else if (j < points_end && i >= points_end && i < images_end)
    {
      Index point_id = j / params_per_point_;
      Index image_id = (i - points_end) / extrinsic_params_per_image_;
      Index row_offset = j % params_per_point_;
      Index col_offset = (i - points_end) % extrinsic_params_per_image_;
      const PointImageBlock* block_ptr = GetPointImageBlock(point_id,
                                                            image_id);
      if (block_ptr == NULL)
      {
        return Scalar(0);
      }
      else
      {
        return (*block_ptr)(row_offset, col_offset);
      }
    }
    else if (i < points_end && j >= images_end && j < cameras_end)
    {
      Index point_id = i / params_per_point_;
      Index camera_id = (j - images_end) / camera_params_size_;
      Index row_offset = i % params_per_point_;
      Index col_offset = (j - images_end) % camera_params_size_;
      const PointCameraBlock* block_ptr = GetPointCameraBlock(point_id,
                                                              camera_id);
      if (block_ptr == NULL)
      {
        return Scalar(0);
      }
      else
      {
        return (*block_ptr)(row_offset, col_offset);
      }
    }
    else if (j < points_end && i >= images_end && i < cameras_end)
    {
      Index point_id = j / params_per_point_;
      Index camera_id = (i - images_end) / camera_params_size_;
      Index row_offset = j % params_per_point_;
      Index col_offset = (i - images_end) % camera_params_size_;
      const PointCameraBlock* block_ptr = GetPointCameraBlock(point_id,
                                                              camera_id);
      if (block_ptr == NULL)
      {
        return Scalar(0);
      }
      else
      {
        return (*block_ptr)(row_offset, col_offset);
      }
    }
    else if (i >= points_end && i < images_end && j >= images_end)
    {
      Index image_id = (i - points_end) / extrinsic_params_per_image_;
      Index camera_id = (j - images_end) / camera_params_size_;
      Index row_offset = (i - points_end) % extrinsic_params_per_image_;
      Index col_offset = (j - images_end) % camera_params_size_;
      const ImageCameraBlock* block_ptr = GetImageCameraBlock(image_id,
                                                              camera_id);
      if (block_ptr == NULL)
      {
        return Scalar(0);
      }
      else
      {
        return (*block_ptr)(row_offset, col_offset);
      }
    }
    else if (j >= points_end && j < images_end && i >= images_end)
    {
      Index image_id = (j - points_end) / extrinsic_params_per_image_;
      Index camera_id = (i - images_end) / camera_params_size_;
      Index row_offset = (j - points_end) % extrinsic_params_per_image_;
      Index col_offset = (i - images_end) % camera_params_size_;
      const ImageCameraBlock* block_ptr = GetImageCameraBlock(image_id,
                                                              camera_id);
      if (block_ptr == NULL)
      {
        return Scalar(0);
      }
      else
      {
        return (*block_ptr)(row_offset, col_offset);
      }
    }
    else
    {
      return Scalar(0);
    }
  }

  void Resize(size_t number_of_points,
              size_t number_of_images,
              size_t number_of_cameras,
              Index camera_params_size)
  {
    number_of_points_ = number_of_points;
    number_of_images_ = number_of_images;
    number_of_cameras_ = number_of_cameras;
    camera_params_size_ = camera_params_size;

    point_blocks_.resize(number_of_points_, PointBlock::Zero());
    image_blocks_.resize(number_of_images_, ImageBlock::Zero());
    camera_blocks_.resize(number_of_cameras_,
                          CameraBlock::Zero(camera_params_size_,
                                            camera_params_size_));

    point_image_map_.resize(Index(number_of_points_),
                            Index(number_of_images_));
    point_camera_map_.resize(Index(number_of_points_),
                             Index(number_of_cameras_));
    image_camera_map_.resize(Index(number_of_images_),
                             Index(number_of_cameras_));

    point_image_blocks_.clear();
    point_camera_blocks_.clear();
    image_camera_blocks_.clear();

  }

  const PointBlock& GetPointBlock(Index point_id) const
  {
    return point_blocks_[point_id];
  }
  PointBlock& GetPointBlock(Index point_id)
  {
    return point_blocks_[point_id];
  }

  const ImageBlock& GetImageBlock(Index image_id) const
  {
    return image_blocks_[image_id];
  }
  ImageBlock& GetImageBlock(Index image_id)
  {
    return image_blocks_[image_id];
  }

  const CameraBlock& GetCameraBlock(Index camera_id) const
  {
    return camera_blocks_[camera_id];
  }
  CameraBlock& GetCameraBlock(Index camera_id)
  {
    return camera_blocks_[camera_id];
  }

  Err AddPointImageBlock(Index point_id,
                          Index image_id,
                          const PointImageBlock& point_image_block)
  {
    BlockIdx block_id = point_image_map_.coeff(point_id, image_id);
    if (block_id == 0)
    {
      return -1;
    }
    else
    {
      point_image_blocks_[block_id - 1] += point_image_block;
      return 0;
    }
  }

  Err AddPointCameraBlock(Index point_id,
                          Index camera_id,
                          const PointCameraBlock& point_camera_block)
  {
    BlockIdx block_id = point_camera_map_.coeff(point_id, camera_id);
    if (block_id == 0)
    {
      return -1;
    }
    else
    {
      point_camera_blocks_[block_id - 1] += point_camera_block;
      return 0;
    }
    return 0;
  }

  Err AddImageCameraBlock(Index image_id,
                           Index camera_id,
                           const ImageCameraBlock& image_camera_block)
  {
    BlockIdx block_id = image_camera_map_.coeff(image_id, camera_id);
    if (block_id == 0)
    {
      return -1;
    }
    else
    {
      image_camera_blocks_[block_id - 1] += image_camera_block;
    }
    return 0;
  }

  template <typename InputIterators>
  void SetPointImageMap(const InputIterators& begin,
                        const InputIterators& end)
  {
    point_image_map_.setFromTriplets(begin, end);
    point_image_blocks_.resize(size_t(point_image_map_.nonZeros()),
                               PointImageBlock::Zero());
  }

  template <typename InputIterators>
  void SetPointCameraMap(const InputIterators& begin,
                         const InputIterators& end)
  {
    point_camera_map_.setFromTriplets(begin, end);
    point_camera_blocks_.resize(size_t(point_camera_map_.nonZeros()),
                             PointCameraBlock::Zero(params_per_point_,
                                                    camera_params_size_));
  }

  template <typename InputIterators>
  void SetImageCameraMap(const InputIterators& begin,
                         const InputIterators& end)
  {
    image_camera_map_.setFromTriplets(begin, end);
    image_camera_blocks_.resize(
      size_t(image_camera_map_.nonZeros()),
      ImageCameraBlock::Zero(extrinsic_params_per_image_,
                             camera_params_size_));
  }

  const PointImageBlock* GetPointImageBlock(Index point_id,
                                            Index image_id) const
  {
    BlockIdx block_id = point_image_map_.coeff(point_id, image_id);
    if (block_id == 0)
    {
      return NULL;
    }
    else
    {
      return &(point_image_blocks_[block_id - 1]);
    }
  }

  const PointCameraBlock* GetPointCameraBlock(Index point_id,
                                              Index camera_id) const
  {
    BlockIdx block_id = point_camera_map_.coeff(point_id, camera_id);
    if (block_id == 0)
    {
      return NULL;
    }
    else
    {
      return &(point_camera_blocks_[block_id - 1]);
    }
  }

  const ImageCameraBlock* GetImageCameraBlock(Index image_id,
                                              Index camera_id) const
  {
    BlockIdx block_id = image_camera_map_.coeff(image_id, camera_id);
    if (block_id == 0)
    {
      return NULL;
    }
    else
    {
      return &(image_camera_blocks_[block_id - 1]);
    }
  }

private:
  size_t number_of_points_;
  size_t number_of_images_;
  size_t number_of_cameras_;
  Index camera_params_size_;
  PointBlockContainer point_blocks_;
  ImageBlockContainer image_blocks_;
  CameraBlockContainer camera_blocks_;
  PointImageBlockContainer point_image_blocks_;
  PointCameraBlockContainer point_camera_blocks_;
  ImageCameraBlockContainer image_camera_blocks_;
  SparseBlockMap point_image_map_;
  SparseBlockMap point_camera_map_;
  SparseBlockMap image_camera_map_;

};

}
}
}

#endif
