#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NORMAL_EQUATION_SOLVER_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_NORMAL_EQUATION_SOLVER_HPP_

#include "hs_sfm/bundle_adjustment/camera_shared_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_augmented_normal_matrix.hpp"
#include "hs_sfm/bundle_adjustment/camera_shared_gradient.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class CameraSharedNormalEquationSolver
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef CameraSharedVectorFunction<Scalar> VectorFunction;
  typedef CameraSharedAugmentedNormalMatrix<Scalar> AugmentedNormalMatrix;
  typedef CameraSharedGradient<Scalar> Gradient;

  typedef typename VectorFunction::XVector DeltaXVector;

private:
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;
  static const Index params_per_point_ = VectorFunction::params_per_point_;
  static const Index params_per_key_ = VectorFunction::params_per_key_;
  static const Index extrinsic_params_per_image_ =
    VectorFunction::extrinsic_params_per_image_;
  typedef typename AugmentedNormalMatrix::NormalMatrix NormalMatrix;
  typedef typename NormalMatrix::PointBlock PointBlock;
  typedef typename NormalMatrix::ImageBlock ImageBlock;
  typedef typename NormalMatrix::CameraBlock CameraBlock;
  typedef typename NormalMatrix::PointImageBlock PointImageBlock;
  typedef typename NormalMatrix::PointCameraBlock PointCameraBlock;
  typedef typename NormalMatrix::ImageCameraBlock ImageCameraBlock;

  typedef typename Gradient::PointSegment PointSegment;

  typedef EIGEN_MATRIX(Scalar, extrinsic_params_per_image_, params_per_point_)
          ImagePointBlock;
  typedef EIGEN_MATRIX(Scalar, Eigen::Dynamic, params_per_point_)
          CameraPointBlock;
  typedef EIGEN_STD_VECTOR(ImagePointBlock) ImagePointBlockContainer;
  typedef EIGEN_STD_VECTOR(CameraPointBlock) CameraPointBlockContainer;
  typedef typename ImagePointBlockContainer::size_type BlockIdx;
  enum
  {
    EigenDefaultMajor =
#if defined(EIGEN_DEFAULT_TO_ROW_MAJOR)
      Eigen::RowMajor
#else
      Eigen::ColMajor
#endif
  };
  typedef EIGEN_SPARSE_MATRIX(BlockIdx, EigenDefaultMajor, Index)
          SparseBlockMap;
  typedef Eigen::Triplet<BlockIdx, Index> BlockTriplet;
  typedef std::vector<BlockTriplet> BlockTripletContainer;

public:
  typedef EIGEN_SPARSE_MATRIX(Scalar, EigenDefaultMajor, Index)
          ReducedMatrix;
  typedef Eigen::Triplet<Scalar, Index> ReducedTriplet;
  typedef std::vector<ReducedTriplet> ReducedTripletContainer;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) ReducedRHS;

public:
  Err operator() (const AugmentedNormalMatrix& augmented_normal_matrix,
                  const Gradient& gradient,
                  DeltaXVector& delta_x) const
  {
    ReducedMatrix reduced_matrix;
    ReducedRHS reduced_rhs;
    if (ComputeReducedSystem(augmented_normal_matrix, gradient,
                             reduced_matrix, reduced_rhs) != 0) return -1;

    const NormalMatrix& normal_matrix =
      augmented_normal_matrix.normal_matrix_ref();
    Scalar mu = augmented_normal_matrix.mu();
    Index number_of_points = Index(normal_matrix.number_of_points());
    Index number_of_images = Index(normal_matrix.number_of_images());
    Index number_of_cameras = Index(normal_matrix.number_of_cameras());
    Index camera_params_size = normal_matrix.camera_params_size();
    Index x_size = number_of_points * params_per_point_ +
                   number_of_images * extrinsic_params_per_image_ +
                   number_of_cameras * camera_params_size;
    delta_x.resize(x_size);
    if (SolveReducedSystem(reduced_matrix, reduced_rhs,
                           number_of_points * params_per_point_,
                           delta_x) != 0)
      return -1;

    if (BackSubstitution(augmented_normal_matrix, gradient, delta_x) != 0)
      return -1;

    return 0;
  }

  Err ComputeReducedSystem(const AugmentedNormalMatrix& augmented_normal_matrix,
                           const Gradient& gradient,
                           ReducedMatrix& reduced_matrix,
                           ReducedRHS& reduced_rhs) const
  {
    const NormalMatrix& normal_matrix =
      augmented_normal_matrix.normal_matrix_ref();
    Scalar mu = augmented_normal_matrix.mu();
    Index number_of_points = Index(normal_matrix.number_of_points());
    Index number_of_images = Index(normal_matrix.number_of_images());
    Index number_of_cameras = Index(normal_matrix.number_of_cameras());
    Index camera_params_size = normal_matrix.camera_params_size();
    Index cameras_begin = number_of_images * extrinsic_params_per_image_;
    Index reduced_size = cameras_begin + number_of_cameras * camera_params_size;

    ImagePointBlockContainer weighted_image_point_blocks;
    BlockTripletContainer image_point_triplets;
    CameraPointBlockContainer weighted_camera_point_blocks;
    BlockTripletContainer camera_point_triplets;
    for (Index point_id = 0; point_id < number_of_points; point_id++)
    {
      PointBlock point_block = normal_matrix.GetPointBlock(point_id);
      for (Index i = 0; i < params_per_point_; i++)
      {
        point_block(i, i) += mu;
      }
      PointBlock point_inverse = point_block.inverse();
      for (Index image_id = 0; image_id < number_of_images; image_id++)
      {
        const PointImageBlock* point_image_block =
          normal_matrix.GetPointImageBlock(point_id, image_id);
        if (point_image_block != NULL)
        {
          ImagePointBlock weighted_image_point_block =
            (*point_image_block).transpose() * point_inverse;
          weighted_image_point_blocks.push_back(weighted_image_point_block);
          image_point_triplets.push_back(
            BlockTriplet(image_id, point_id,
                         weighted_image_point_blocks.size()));
        }
      }//for (Index image_id = 0; image_id < number_of_images; image_id++)

      for (Index camera_id = 0; camera_id < number_of_cameras; camera_id++)
      {
        const PointCameraBlock* point_camera_block =
          normal_matrix.GetPointCameraBlock(point_id, camera_id);
        if (point_camera_block != NULL)
        {
          CameraPointBlock weighted_camera_point_block =
            (*point_camera_block).transpose() * point_inverse;
          weighted_camera_point_blocks.push_back(weighted_camera_point_block);
          camera_point_triplets.push_back(
            BlockTriplet(camera_id, point_id,
                         weighted_camera_point_blocks.size()));
        }
      }//for (Index camera_id = 0; camera_id < number_of_cameras; camera_id++)
    }//for (Index point_id = 0; point_id < number_of_points; point_id++
    SparseBlockMap image_point_map(number_of_images, number_of_points);
    image_point_map.setFromTriplets(image_point_triplets.begin(),
                                    image_point_triplets.end());
    SparseBlockMap camera_point_map(number_of_cameras, number_of_points);
    camera_point_map.setFromTriplets(camera_point_triplets.begin(),
                                     camera_point_triplets.end());

    ReducedTripletContainer reduced_triplets;
    for (Index weighted_image_id = 0; weighted_image_id < number_of_images;
         weighted_image_id++)
    {
      for (Index image_id = weighted_image_id;
           image_id < number_of_images; image_id++)
      {
        ImageBlock image_block = ImageBlock::Zero();
        if (weighted_image_id == image_id)
        {
          image_block = normal_matrix.GetImageBlock(image_id);
          for (Index i = 0; i < extrinsic_params_per_image_; i++)
          {
            image_block(i, i) += mu;
          }
        }
        bool exist = false;
        for (Index point_id = 0; point_id < number_of_points; point_id++)
        {
          BlockIdx weighted_image_point_id =
            image_point_map.coeff(weighted_image_id, point_id);
          const PointImageBlock* point_image_block =
            normal_matrix.GetPointImageBlock(point_id, image_id);
          if (weighted_image_point_id > 0 && point_image_block != NULL)
          {
            const ImagePointBlock& weighted_image_point_block =
              weighted_image_point_blocks[weighted_image_point_id - 1];
            image_block -= weighted_image_point_block * (*point_image_block);
            exist = true;
          }
        }// for (Index point_id = 0; point_id < number_of_points; point_id++)

        if (exist)
        {
          if (weighted_image_id == image_id)
          {
            for (Index i = 0; i < extrinsic_params_per_image_; i++)
            {
              for (Index j = 0; j < extrinsic_params_per_image_; j++)
              {
                reduced_triplets.push_back(
                  ReducedTriplet(image_id * extrinsic_params_per_image_ + i,
                                 image_id * extrinsic_params_per_image_ + j,
                                 image_block(i, j)));
              }// for (Index j = 0; j < extrinsic_params_per_image_; j++)
            }// for (Index i = 0; i < extrinsic_params_per_image_; i++)
          }// if (weighted_image_id == image_id)
          else
          {
            for (Index i = 0; i < extrinsic_params_per_image_; i++)
            {
              for (Index j = 0; j < extrinsic_params_per_image_; j++)
              {
                Index row_id =
                  weighted_image_id * extrinsic_params_per_image_ + i;
                Index col_id =
                  image_id * extrinsic_params_per_image_ + j;
                reduced_triplets.push_back(
                  ReducedTriplet(row_id, col_id, image_block(i, j)));
                reduced_triplets.push_back(
                  ReducedTriplet(col_id, row_id, image_block(i, j)));
              }// for (Index j = 0; j < extrinsic_params_per_image_; j++)
            }// for (Index i = 0; i < extrinsic_params_per_image_; i++)
          }
        }// if (exist)
      }
    }

    for (Index weighted_camera_id = 0;
         weighted_camera_id < number_of_cameras; weighted_camera_id++)
    {
      for (Index camera_id = weighted_camera_id;
           camera_id < number_of_cameras; camera_id++)
      {
        CameraBlock camera_block = CameraBlock::Zero(camera_params_size,
                                                     camera_params_size);
        if (weighted_camera_id == camera_id)
        {
          camera_block = normal_matrix.GetCameraBlock(camera_id);
          for (Index i = 0; i < camera_params_size; i++)
          {
            camera_block(i, i) += mu;
          }
        }
        bool exist = false;
        for (Index point_id = 0; point_id < number_of_points; point_id++)
        {
          BlockIdx weighted_camera_point_id =
            camera_point_map.coeff(weighted_camera_id, point_id);
          const PointCameraBlock* point_camera_block =
            normal_matrix.GetPointCameraBlock(point_id, camera_id);
          if (weighted_camera_point_id > 0 && point_camera_block != NULL)
          {
            const CameraPointBlock& weighted_camera_point_block =
              weighted_camera_point_blocks[weighted_camera_point_id - 1];
            camera_block -= weighted_camera_point_block * (*point_camera_block);
            exist = true;
          }
        }//for (Index point_id = 0; point_id < number_of_points; point_id++)

        if (exist)
        {
          if (weighted_camera_id == camera_id)
          {
            for (Index i = 0; i < camera_params_size; i++)
            {
              for (Index j = 0; j < camera_params_size; j++)
              {
                Index row_id =
                  cameras_begin + camera_id * camera_params_size + i;
                Index col_id =
                  cameras_begin + camera_id * camera_params_size + j;
                reduced_triplets.push_back(
                  ReducedTriplet(row_id, col_id, camera_block(i, j)));
              }//for (Index j = 0; j < camera_params_size; j++)
            }// for (Index i = 0; i < camera_params_size; i++)
          }// if (weighted_camera_id == camera_id)
          else
          {
            for (Index i = 0; i < camera_params_size; i++)
            {
              for (Index j = 0; j < camera_params_size; j++)
              {
                Index row_id =
                  cameras_begin + weighted_camera_id * camera_params_size + i;
                Index col_id =
                  cameras_begin + camera_id * camera_params_size + j;
                reduced_triplets.push_back(
                  ReducedTriplet(row_id, col_id, camera_block(i, j)));
                reduced_triplets.push_back(
                  ReducedTriplet(col_id, row_id, camera_block(i, j)));
              }//for (Index j = 0; j < camera_params_size; j++)
            }// for (Index i = 0; i < camera_params_size; i++)
          }
        }// if (exist)
      }
    }

    for (Index image_id = 0; image_id < number_of_images; image_id++)
    {
      for (Index camera_id = 0; camera_id < number_of_cameras; camera_id++)
      {
        ImageCameraBlock image_camera_block =
          ImageCameraBlock::Zero(extrinsic_params_per_image_,
                                 camera_params_size);
        const ImageCameraBlock* image_camera_block_ptr =
          normal_matrix.GetImageCameraBlock(image_id, camera_id);
        if (image_camera_block_ptr != NULL)
        {
          image_camera_block = *image_camera_block_ptr;
        }
        bool exist = false;
        for (Index point_id = 0; point_id < number_of_points; point_id++)
        {
          const PointCameraBlock* point_camera_block =
            normal_matrix.GetPointCameraBlock(point_id, camera_id);
          BlockIdx weighted_image_point_id =
            image_point_map.coeff(image_id, point_id);
          if (weighted_image_point_id > 0 && point_camera_block != NULL)
          {
            const ImagePointBlock& weighted_image_point_block =
              weighted_image_point_blocks[weighted_image_point_id - 1];
            image_camera_block -= weighted_image_point_block *
                                  (*point_camera_block);
            exist = true;
          }// if (weighted_point_camera_id > 0 && point_image_block != NULL)
        }// for (Index point_id = 0; point_id < number_of_points; point_id++)
        if (exist)
        {
          for (Index i = 0; i < extrinsic_params_per_image_; i++)
          {
            for (Index j = 0; j < camera_params_size; j++)
            {
              Index row_id = image_id * extrinsic_params_per_image_ + i;
              Index col_id = cameras_begin +
                             camera_id * camera_params_size + j;
              reduced_triplets.push_back(
                ReducedTriplet(row_id, col_id, image_camera_block(i, j)));
              reduced_triplets.push_back(
                ReducedTriplet(col_id, row_id, image_camera_block(i, j)));
            }// for (Index j = 0; j < camera_params_size; j++)
          }// for (Index i = 0; i < extrinsic_params_per_image_; i++)
        }
      }// for (Index camera_id = 0; camera_id < number_of_cameras; camera_id++)
    }// for (Index image_id = 0; image_id < number_of_images; image_id++)

    reduced_matrix.resize(reduced_size, reduced_size);
    reduced_matrix.setFromTriplets(reduced_triplets.begin(),
                                   reduced_triplets.end());

    reduced_rhs.resize(reduced_size);
    for (Index image_id = 0; image_id < number_of_images; image_id++)
    {
      reduced_rhs.segment(image_id * extrinsic_params_per_image_,
                          extrinsic_params_per_image_) =
        gradient.GetImageSegment(image_id);
      for (Index point_id = 0; point_id < number_of_points; point_id++)
      {
        BlockIdx weighted_image_point_id =
          image_point_map.coeff(image_id, point_id);
        if (weighted_image_point_id > 0)
        {
          const ImagePointBlock& weighted_image_point_block =
            weighted_image_point_blocks[weighted_image_point_id - 1];
          reduced_rhs.segment(image_id * extrinsic_params_per_image_,
                              extrinsic_params_per_image_) -=
            weighted_image_point_block * gradient.GetPointSegment(point_id);
        }
      }// for (Index point_id = 0; point_id < number_of_points; point_id++)
    }// for (Index image_id = 0; image_id < number_of_images; image_id++)

    for (Index camera_id = 0; camera_id < number_of_cameras; camera_id++)
    {
      reduced_rhs.segment(cameras_begin + camera_id * camera_params_size,
                          camera_params_size) =
        gradient.GetCameraSegment(camera_id);
      for (Index point_id; point_id < number_of_points; point_id++)
      {
        BlockIdx weighted_camera_point_id =
          camera_point_map.coeff(camera_id, point_id);
        if (weighted_camera_point_id > 0)
        {
          const CameraPointBlock& weighted_camera_point_block =
            weighted_camera_point_blocks[weighted_camera_point_id - 1];
          reduced_rhs.segment(cameras_begin + camera_id * camera_params_size,
                              camera_params_size) -=
            weighted_camera_point_block * gradient.GetPointSegment(point_id);
        }
      }// for (Index point_id; point_id < number_of_points; point_id++)
    }// for (Index camera_id = 0; camera_id < number_of_cameras; camera_id++)

    return 0;
  }

  Err SolveReducedSystem(const ReducedMatrix& reduced_matrix,
                         const ReducedRHS& reduced_rhs,
                         Index reduced_system_begin,
                         DeltaXVector& delta_x) const
  {
    Eigen::SimplicialLDLT<ReducedMatrix> solver;
    solver.compute(reduced_matrix);
    if (solver.info() != Eigen::Success)
    {
      return -1;
    }
    delta_x.segment(reduced_system_begin, reduced_rhs.rows()) =
      solver.solve(reduced_rhs);
    if (solver.info() != Eigen::Success)
    {
      return -1;
    }
    return 0;
  }

  Err BackSubstitution(const AugmentedNormalMatrix& augmented_normal_matrix,
                       const Gradient& gradient,
                       DeltaXVector& delta_x) const
  {
    const NormalMatrix& normal_matrix =
      augmented_normal_matrix.normal_matrix_ref();
    Scalar mu = augmented_normal_matrix.mu();
    Index number_of_points = Index(normal_matrix.number_of_points());
    Index number_of_images = Index(normal_matrix.number_of_images());
    Index number_of_cameras = Index(normal_matrix.number_of_cameras());
    Index camera_params_size = normal_matrix.camera_params_size();
    Index cameras_begin = number_of_images * extrinsic_params_per_image_;
    Index reduced_size = cameras_begin + number_of_cameras * camera_params_size;

    for (Index point_id = 0; point_id < number_of_points; point_id++)
    {
      PointBlock point_block = normal_matrix.GetPointBlock(point_id);
      for (Index i = 0; i < params_per_point_; i++)
      {
        point_block(i, i) += mu;
      }
      PointBlock point_block_inverse = point_block.inverse();
      PointSegment point_segment = point_block_inverse *
                                   gradient.GetPointSegment(point_id);
      for (Index image_id = 0; image_id < number_of_images; image_id++)
      {
        const PointImageBlock* point_image_block =
          normal_matrix.GetPointImageBlock(point_id, image_id);
        if (point_image_block != NULL)
        {
          point_segment -=
            point_block_inverse * (*point_image_block) *
            delta_x.segment(number_of_points * params_per_point_ +
                            image_id * extrinsic_params_per_image_,
                            extrinsic_params_per_image_);
        }
      }

      for (Index camera_id = 0; camera_id < number_of_cameras; camera_id++)
      {
        const PointCameraBlock* point_camera_block =
          normal_matrix.GetPointCameraBlock(point_id, camera_id);
        if (point_camera_block != NULL)
        {
          point_segment -=
            point_block_inverse * (*point_camera_block) *
            delta_x.segment(number_of_points * params_per_point_ +
                            cameras_begin +
                            camera_id * camera_params_size,
                            camera_params_size);
        }
      }

      delta_x.segment(point_id * params_per_point_, params_per_point_) =
        point_segment;
    }

    return 0;
  }

};

}
}
}

#endif
