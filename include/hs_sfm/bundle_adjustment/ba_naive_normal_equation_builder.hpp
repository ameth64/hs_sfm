#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_EQUATION_BUILDER_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_NORMAL_EQUATION_BUILDER_HPP_

#include <cassert>
#include <cmath>

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_feature_covariance_inverse.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_jacobian_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_normal_matrix.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_gradient.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BANaiveNormalEquationBuilder
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef BANaiveVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::Index Index;
  typedef BANaiveJacobianMatrix<Scalar, Index,
                                VectorFunction::params_per_feature_,
                                VectorFunction::params_per_camera_,
                                VectorFunction::params_per_point_>
          JacobianMatrix;
  typedef BANaiveFeatureCovarianceInverse<Scalar> YCovarianceInverse;

  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) Residuals;

  typedef BANaiveNormalMatrix<Scalar, Index,
                              VectorFunction::params_per_camera_,
                              VectorFunction::params_per_point_> NormalMatrix;
  typedef BANaiveGradient<Scalar, Index,
                          VectorFunction::params_per_camera_,
                          VectorFunction::params_per_point_> Gradient;

  typedef typename JacobianMatrix::CameraDerivativeBlock CamDrvBlk;
  typedef typename JacobianMatrix::PointDerivativeBlock PtDrvBlk;
  typedef typename JacobianMatrix::DerivativeId DerivativeId;
  typedef typename JacobianMatrix::DerivativeMap DrivativeMap;
  typedef typename YCovarianceInverse::Block YCovInvBlk;
  typedef typename NormalMatrix::CameraBlock CameraBlock;
  typedef typename NormalMatrix::PointBlock PointBlock;
  typedef typename NormalMatrix::MixBlock MixBlock;
  typedef typename NormalMatrix::MixBlockIdx MixBlockIdx;
  typedef typename NormalMatrix::MixBlockMap MixBlockMap;
  typedef typename Gradient::CameraSegment CameraSegment;
  typedef typename Gradient::PointSegment PointSegment;

  Err operator()(const JacobianMatrix& jacobian_matrix,
                 const Residuals& residuals,
                 const YCovarianceInverse& y_covariance_inverse,
                 NormalMatrix& normal_matrix,
                 Gradient& gradient) const
  {
    if (BuildNormalMatrix(jacobian_matrix,
                          y_covariance_inverse,
                          normal_matrix) != 0) return -1;
    if (BuildGradient(jacobian_matrix,
                     residuals,
                     y_covariance_inverse,
                     gradient) != 0) return -1;

    return 0;
  }

  Err BuildNormalMatrix(const JacobianMatrix& jacobian_matrix,
                        const YCovarianceInverse& y_covariance_inverse,
                        NormalMatrix& normal_matrix) const
  {
    //Camera blocks
    Index number_of_cameras = jacobian_matrix.number_of_cameras();
    Index number_of_points = jacobian_matrix.number_of_points();
    normal_matrix.camera_blocks.resize(number_of_cameras, CameraBlock::Zero());
    normal_matrix.point_blocks.resize(number_of_points, PointBlock::Zero());
    size_t number_of_camera_blocks =
      jacobian_matrix.camera_derivatives().size();
    for (size_t i = 0; i < number_of_camera_blocks; i++)
    {
      const CamDrvBlk& cam_drv_blk = jacobian_matrix.camera_derivatives()[i];
      Index feature_id = cam_drv_blk.feature_id;
      const YCovInvBlk& y_cov_inv_blk =
        y_covariance_inverse.blocks[feature_id];
      normal_matrix.camera_blocks[cam_drv_blk.camera_id] +=
        cam_drv_blk.derivative_block.transpose() * y_cov_inv_blk *
        cam_drv_blk.derivative_block;
    }

    size_t number_of_point_blocks =
      jacobian_matrix.point_derivatives().size();
    for (size_t i = 0; i < number_of_point_blocks; i++)
    {
      const PtDrvBlk& pt_drv_blk = jacobian_matrix.point_derivatives()[i];
      Index feature_id = pt_drv_blk.feature_id;
      const YCovInvBlk& y_cov_inv_blk =
        y_covariance_inverse.blocks[feature_id];
      normal_matrix.point_blocks[pt_drv_blk.point_id] +=
        pt_drv_blk.derivative_block.transpose() * y_cov_inv_blk *
        pt_drv_blk.derivative_block;
    }

    typedef Eigen::Triplet<MixBlockIdx, Index> TripletType;
    std::vector<TripletType> mix_block_triplets;
    normal_matrix.mix_blocks.clear();
    normal_matrix.number_of_cameras = number_of_cameras;
    normal_matrix.number_of_points = number_of_points;
    for (Index i = 0; i < number_of_cameras; i++)
    {
      for (Index k = 0; k < number_of_points; k++)
      {
        DerivativeId camera_id =
          jacobian_matrix.camera_derivatives_map().coeff(i, k);
        DerivativeId point_id =
          jacobian_matrix.point_derivatives_map().coeff(i, k);
        if (camera_id != 0 &&
            point_id != 0)
        {
          const CamDrvBlk& cam_drv_blk =
            jacobian_matrix.camera_derivatives()[camera_id - 1];
          const PtDrvBlk& pt_drv_blk =
            jacobian_matrix.point_derivatives()[point_id - 1];
          assert(cam_drv_blk.feature_id == pt_drv_blk.feature_id);
          Index feature_id = cam_drv_blk.feature_id;
          const YCovInvBlk& y_cov_inv_blk =
            y_covariance_inverse.blocks[feature_id];
          MixBlock mix_block =
            cam_drv_blk.derivative_block.transpose() *
            y_cov_inv_blk *
            pt_drv_blk.derivative_block;
          normal_matrix.mix_blocks.push_back(mix_block);
          mix_block_triplets.push_back(
            TripletType(i, k, normal_matrix.mix_blocks.size()));
        }
      }
    }
    normal_matrix.mix_block_map.resize(number_of_cameras, number_of_points);
    normal_matrix.mix_block_map.setFromTriplets(
      mix_block_triplets.begin(), mix_block_triplets.end());

    return 0;
  }

  Err BuildGradient(const JacobianMatrix& jacobian_matrix,
                   const Residuals& residuals,
                   const YCovarianceInverse& y_covariance_inverse,
                   Gradient& gradient) const
  {
    //Camera blocks
    Index number_of_cameras = jacobian_matrix.number_of_cameras();
    Index number_of_points = jacobian_matrix.number_of_points();
    gradient.camera_segments.resize(number_of_cameras, CameraSegment::Zero());
    gradient.point_segments.resize(number_of_points, PointSegment::Zero());
    size_t number_of_camera_blocks =
      jacobian_matrix.camera_derivatives().size();
    for (size_t i = 0; i < number_of_camera_blocks; i++)
    {
      const CamDrvBlk& cam_drv_blk = jacobian_matrix.camera_derivatives()[i];
      Index feature_id = cam_drv_blk.feature_id;
      const YCovInvBlk& y_cov_inv_blk =
        y_covariance_inverse.blocks[feature_id];
      gradient.camera_segments[cam_drv_blk.camera_id] +=
        cam_drv_blk.derivative_block.transpose() * y_cov_inv_blk *
        residuals.segment(VectorFunction::params_per_feature_ * feature_id,
                          VectorFunction::params_per_feature_);
    }

    size_t number_of_point_blocks =
      jacobian_matrix.point_derivatives().size();
    for (size_t i = 0; i < number_of_point_blocks; i++)
    {
      const PtDrvBlk& pt_drv_blk = jacobian_matrix.point_derivatives()[i];
      Index feature_id = pt_drv_blk.feature_id;
      const YCovInvBlk& y_cov_inv_blk =
        y_covariance_inverse.blocks[feature_id];
      gradient.point_segments[pt_drv_blk.point_id] +=
        pt_drv_blk.derivative_block.transpose() * y_cov_inv_blk *
        residuals.segment(VectorFunction::params_per_feature_ * feature_id,
                          VectorFunction::params_per_feature_);
    }

    gradient.number_of_cameras = number_of_cameras;
    gradient.number_of_points = number_of_points;

    return 0;
  }
};

}
}
}

#endif
