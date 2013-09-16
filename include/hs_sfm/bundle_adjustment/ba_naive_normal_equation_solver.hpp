#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_EQUATION_SOLVER_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_EQUATION_SOLVER_HPP_

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_augmented_normal_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_gradient.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BANaiveNormalEquationSolver
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;
  typedef BANaiveAugmentedNormalMatrix<Scalar, Index,
                                       VectorFunction::params_per_camera_,
                                       VectorFunction::params_per_point_> 
          AugmentedNormalMatrix;
  typedef typename AugmentedNormalMatrix::NormalMatrix NormalMatrix;
  typedef BANaiveGradient<Scalar, Index,
                          VectorFunction::params_per_camera_,
                          VectorFunction::params_per_point_> Gradient;

  typedef typename NormalMatrix::CameraBlock CameraBlock;
  typedef typename NormalMatrix::PointBlock PointBlock;
  typedef typename NormalMatrix::MixBlock MixBlock;
  typedef typename NormalMatrix::MixBlockIdx MixBlockIdx;
  typedef typename NormalMatrix::MixBlockMap MixBlockMap;
  typedef typename Gradient::CameraSegment CameraSegment;
  typedef typename Gradient::PointSegment PointSegment;

  enum
  {
    EigenDefaultMajor =
#if defined(EIGEN_DEFAULT_TO_ROW_MAJOR)
      Eigen::RowMajor
#else
      Eigen::ColMajor
#endif
  };
  typedef EIGEN_SPARSE_MATRIX(Scalar, EigenDefaultMajor, Index)
          SchurMatrix;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) CameraParams;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) PointParams;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) CameraRHS;
  typedef EIGEN_MATRIX(Scalar,
                    VectorFunction::params_per_camera_,
                    VectorFunction::params_per_point_) MixPointBlock;
  typedef EIGEN_STD_VECTOR(MixPointBlock) MixPointBlockContainer;
  typedef typename MixPointBlockContainer::size_type MixPointBlockIdx;
  typedef EIGEN_SPARSE_MATRIX(MixPointBlockIdx, EigenDefaultMajor, Index)
          MixPointBlockMap;
  typedef XVector DeltaXVector;

  Err operator()(const AugmentedNormalMatrix& augmented_normal_matrix,
                 const Gradient& gradient,
                 DeltaXVector& delta_x) const
  {
    if (SolveCameraParams(augmented_normal_matrix, gradient, delta_x) != 0) return -1;
    if (SolvePointParams(augmented_normal_matrix, gradient, delta_x) != 0) return -1;
    return 0;
  }

private:
  Err SolveCameraParams(const AugmentedNormalMatrix& augmented_normal_matrix,
                        const Gradient& gradient,
                        DeltaXVector& delta_x) const
  {
    const NormalMatrix& normal_matrix =
      augmented_normal_matrix.normal_matrix_ref_;
    Index schur_size = normal_matrix.number_of_cameras *
                       VectorFunction::params_per_camera_;
    SchurMatrix schur_mat;
    delta_x.resize(schur_size + normal_matrix.number_of_points *
                                VectorFunction::params_per_point_);

    //compute mix point blocks
    MixPointBlockContainer mix_point_blocks;
    MixPointBlockMap mix_point_block_map;
    typedef Eigen::Triplet<MixPointBlockIdx, Index> MixPointTripletType;
    std::vector<MixPointTripletType> mix_point_block_triplets;
    for (Index i = 0; i < normal_matrix.number_of_cameras; i++)
    {
      for (Index j = 0; j < normal_matrix.number_of_points; j++)
      {
        MixBlockIdx mix_id =
          normal_matrix.mix_block_map.coeff(i, j);
        if (mix_id > 0)
        {
          const MixBlock& mix =
            normal_matrix.mix_blocks[mix_id - 1];
          PointBlock point =
            normal_matrix.point_blocks[j];
          for (Index k = 0; k < VectorFunction::params_per_point_; k++)
          {
            point(k, k) += augmented_normal_matrix.mu_;
          }
          MixPointBlock mix_point = mix * point.inverse();
          mix_point_blocks.push_back(mix_point);
          mix_point_block_triplets.push_back(
            MixPointTripletType(i, j, mix_point_blocks.size()));
        }
      }
    }
    mix_point_block_map.resize(normal_matrix.number_of_cameras, 
                               normal_matrix.number_of_points);
    mix_point_block_map.setFromTriplets(mix_point_block_triplets.begin(),
                                        mix_point_block_triplets.end());

    //compute schur matrix
    typedef Eigen::Triplet<Scalar, Index> SchurTripletType;
    std::vector<SchurTripletType> schur_triplets;
    for (Index i = 0; i < normal_matrix.number_of_cameras; i++)
    {
      for (Index j = i; j < normal_matrix.number_of_cameras; j++)
      {
        CameraBlock camera_block = CameraBlock::Zero();
        for (Index k = 0; k < normal_matrix.number_of_points; k++)
        {
          MixPointBlockIdx mix_point_id = mix_point_block_map.coeff(i, k);
          MixBlockIdx mix_id = normal_matrix.mix_block_map.coeff(j, k);
          if (mix_point_id != 0 &&  mix_id != 0)
          {
            camera_block -=
              mix_point_blocks[mix_point_id - 1] *
              normal_matrix.mix_blocks[mix_id - 1].transpose();
          }
        }

        //fill schur matrix
        if (i == j)
        {
          camera_block += normal_matrix.camera_blocks[i];
          for (Index m = 0; m < VectorFunction::params_per_camera_; m++)
          {
            camera_block(m, m) += augmented_normal_matrix.mu_;
          }

          for (Index m = 0; m < VectorFunction::params_per_camera_; m++)
          {
            for (Index n = 0; n < VectorFunction::params_per_camera_; n++)
            {
              schur_triplets.push_back(
                SchurTripletType(i * VectorFunction::params_per_camera_ + m,
                                 j * VectorFunction::params_per_camera_ + n,
                                 camera_block(m, n)));
            }
          }
        }
        else
        {
          for (Index m = 0; m < VectorFunction::params_per_camera_; m++)
          {
            for (Index n = 0; n < VectorFunction::params_per_camera_; n++)
            {
              schur_triplets.push_back(
                SchurTripletType(i * VectorFunction::params_per_camera_ + m,
                                 j * VectorFunction::params_per_camera_ + n,
                                 camera_block(m, n)));
              schur_triplets.push_back(
                SchurTripletType(j * VectorFunction::params_per_camera_ + n,
                                 i * VectorFunction::params_per_camera_ + m,
                                 camera_block(m, n)));
            }
          }
        }
      }
    }
    schur_mat.resize(schur_size, schur_size);
    schur_mat.setFromTriplets(schur_triplets.begin(),
                              schur_triplets.end());

    //compute camera RHS
    CameraRHS camera_rhs(schur_size);
    for (Index i = 0; i < normal_matrix.number_of_cameras; i++)
    {
      CameraSegment camera_segment = gradient.camera_segments[i];
      for (Index j = 0; j < normal_matrix.number_of_points; j++)
      {
        MixPointBlockIdx mix_point_id = mix_point_block_map.coeff(i, j);
        if (mix_point_id > 0)
        {
          const PointSegment& point_segment = gradient.point_segments[j];
          camera_segment -= mix_point_blocks[mix_point_id - 1] *
                         point_segment;
        }
      }
      camera_rhs.segment(i * VectorFunction::params_per_camera_,
                         VectorFunction::params_per_camera_) = camera_segment;
    }

    //solve camera params
    typedef Eigen::SimplicialLDLT<SchurMatrix> Solver;
    Solver solver;
    solver.compute(schur_mat);
    if (solver.info() != Eigen::Success)
    {
      return -1;
    }
    delta_x.segment(0, schur_size) = solver.solve(camera_rhs);
    if (solver.info() != Eigen::Success)
    {
      return -1;
    }

    return 0;
  }

  Err SolvePointParams(const AugmentedNormalMatrix& augmented_normal_matrix, 
                       const Gradient& gradient,
                       DeltaXVector& delta_x) const
  {
    const NormalMatrix& normal_matrix =
      augmented_normal_matrix.normal_matrix_ref_;
    Index camera_params_size = normal_matrix.number_of_cameras * 
                               VectorFunction::params_per_camera_;
    for (Index i = 0; i < normal_matrix.number_of_points; i++)
    {
      PointSegment point_segment = gradient.point_segments[i];
      for (Index j = 0; j < normal_matrix.number_of_cameras; j++)
      {
        MixBlockIdx mix_id = normal_matrix.mix_block_map.coeff(j, i);
        if (mix_id > 0)
        {
          const MixBlock& mix = normal_matrix.mix_blocks[mix_id - 1];
          point_segment -=
            mix.transpose() *
            delta_x.segment(j * VectorFunction::params_per_camera_,
                            VectorFunction::params_per_camera_);
        }
      }
      PointBlock point_block = normal_matrix.point_blocks[i];
      for (Index j = 0; j < VectorFunction::params_per_point_; j++)
      {
        point_block(j, j) += augmented_normal_matrix.mu_;
      }
      delta_x.segment(camera_params_size +
                      i * VectorFunction::params_per_point_,
                      VectorFunction::params_per_point_) =
        point_block.inverse() * point_segment;
    }

    return 0;
  }
};

}
}
}

#endif
