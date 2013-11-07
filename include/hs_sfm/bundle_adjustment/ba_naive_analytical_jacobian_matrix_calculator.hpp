#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_ANALYTICAL_JACOBIAN_MATRIX_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_NAIVE_ANALYTICAL_JACOBIAN_MATRIX_CALCULATOR_HPP_

#include <vector>

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"
#include "hs_sfm/bundle_adjustment/ba_naive_jacobian_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VectorFunction>
class BANaiveAnalyticalJacobianMatrixCalculator;

/**
 *  计算简单Bundle Adjustment函数的Jacbobian矩阵。
 *
 *  设\f$\mathbf{y_i}\f$为第i个影像特征点对应\f$\mathbf{c_j}\f$为第j个相机的参数，
 *  \f$\mathbf{p_k}\f$为第k个三维点，
 *  则\f$\mathbf{y_i}\f$的齐次表示\f$\mathbf{\hat{y}_i}\f$有：
 *  \f[
 *    \mathbf{\hat{y}_i} = \mathbf{r_j}(\mathbf{p_k}) + \mathbf{t_j}
 *  \f]
 *  其中\f$\mathbf{t_j}\f$表示第j个相机的位移。
 *  \f$\mathbf{r_j}(\mathbf{x})\f$表示第j个相机的朝向，使用轴角旋转表示为：
 *  \f[
 *    \mathbf{r_j}(\mathbf{x}) = cos(\theta)\mathbf{x} +
 *      sin(\theta)((\mathbf{r_j} / \theta) \times \mathbf{x}) +
 *      \frac{1-cos(\theta)}{\theta^2}(\mathbf{r_j}\cdot\mathbf{x})\mathbf{r_j}
 *  \f]
 *  其中\f$\theta\f$为\f$\mathbf{r_j}\f$的模。
 *  于是有\f$\mathbf{\hat{y}_i} = [\hat{y}_{i0},\hat{y}_{i1},\hat{y}_{i2}]^T\f$：
 *  \f[
 *    \hat{y}_{i0} = cos(\theta)p_{k0} + 
 *           \frac{sin(\theta)}{\theta}r_{j1}p_{k2} +
 *           d\frac{1-cos(\theta)}{\theta^2}r_{j0} + t_{j0}
 *  \f]
 *  \f[
 *    \hat{y}_{i1} = cos(\theta)p_{k1} + 
 *           \frac{sin(\theta)}{\theta}r_{j2}p_{k0} +
 *           d\frac{1-cos(\theta)}{\theta^2}r_{j1} + t_{j1}
 *  \f]
 *  \f[
 *    \hat{y}_{i2} = cos(\theta)p_{k2} + 
 *           \frac{sin(\theta)}{\theta}r_{j0}p_{k1} +
 *           d\frac{1-cos(\theta)}{\theta^2}r_{j2} + t_{j2}
 *  \f]
 *  其中\f$d=\mathbf{r_j}\cdot\mathbf{p_k}\f$
 */
template <typename _Scalar>
class BANaiveAnalyticalJacobianMatrixCalculator<BANaiveVectorFunction<_Scalar> >
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

  Err operator()(const VectorFunction& vector_function, const XVector& x,
                 JacobianMatrix& jacobian_matrix) const
  {
    typedef EIGEN_VECTOR(Scalar, 3) Vector3;
    typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;
    typedef Eigen::Triplet<DerivativeId, Index> TripletType;
    Index number_of_cameras = vector_function.number_of_cameras();
    Index number_of_points = vector_function.number_of_points();
    Index number_of_features = vector_function.number_of_features();
    Index camera_params_size = vector_function.GetCameraParamsSize();
    const FeatureMapContainer& feature_maps = vector_function.feature_maps();
    Index x_size = vector_function.GetXSize();
    std::vector<TripletType> cameras_triplets;
    std::vector<TripletType> points_triplets;
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

      Scalar theta = r.norm();
      if (theta == Scalar(0))
      {
        return -1;
      }
      Scalar theta2 = theta * theta;
      Scalar theta4 = theta2 * theta2;
      Scalar d = r.dot(p);

      //sin(theta)
      Scalar st = sin(theta);
      //cos(theta)
      Scalar ct = cos(theta);

      //f = cos(theta) * p
      Vector3 f = ct * p;
      //\frac{\partial f}{\partial r}
      Matrix33 pfpr = p * (-st / theta * r).transpose();
      //\frac{\partial f}{\partial p}
      Matrix33 pfpp = Matrix33::Identity() * ct;

      //s = sin(theta) / theta
      Scalar s = st / theta;
      //\frac{\partial s}{\partial r}
      Vector3 pspr = (ct * theta - st) / theta2 / theta * r;

      //v = r x p
      Vector3 v = r.cross(p);
      //\frac{\partial v}{\partial r}
      Matrix33 pvpr;
      pvpr << 0, p[2], -p[1],
          -p[2], 0, p[0],
          p[1], -p[0], 0;
      //\frac{\partial v}{\partial p}
      Matrix33 pvpp;
      pvpp << 0, -r[2], r[1],
          r[2], 0, -r[0],
          -r[1], r[0], 0;

      //c = (1 - cos(theta)) / theta^2
      Scalar c = (1 - ct) / theta2;
      //\frac{partial c}{\partial r}
      Vector3 pcpr = 
        (st * theta - 2 * (1 - ct)) / theta4 * r;

      //u = r^T * p * r
      Vector3 u = d * r;
      //\frac{partial u}{\partial r}
      Matrix33 pupr;
      pupr << d + r[0] * p[0], r[0] * p[1], r[0] * p[2],
          r[1] * p[0], d + r[1] * p[1], r[1] * p[2],
          r[2] * p[0], r[2] * p[1], d + r[2] * p[2];
      //\frac{partial u}{\patial p}
      Matrix33 pupp;
      pupp << r[0] * r[0], r[0] * r[1], r[0] * r[2],
          r[1] * r[0], r[1] * r[1], r[1] * r[2],
          r[2] * r[0], r[2] * r[1], r[2] * r[2];

      Vector3 y = f + s * v + c * u + t;
      //\frac{partial y}{\partial r}
      Matrix33 pypr = pfpr +
             v * pspr.transpose() + s * pvpr +
             u * pcpr.transpose() + c * pupr;
      //\frac{partial y}{\partial p}
      Matrix33 pypp = pfpp + s * pvpp + c * pupp;
      //\frac{partial y}{\partial t} is identity
      Matrix33 pypt = Matrix33::Identity();

      typename JacobianMatrix::CameraDerivativeBlock camera_block;
      camera_block.camera_id = j;
      camera_block.feature_id = i;
      for (Index m = 0; m < 3; m++)
      {
        //for r
        camera_block.derivative_block(0, m) =
          (pypr(0, m) * y[2] - y[0] * pypr(2, m)) / y[2] / y[2];
        camera_block.derivative_block(1, m) =
          (pypr(1, m) * y[2] - y[1] * pypr(2, m)) / y[2] / y[2];

        //for t
        camera_block.derivative_block(0, 3 + m) =
          (pypt(0, m) * y[2] - y[0] * pypt(2, m)) / y[2] / y[2];
        camera_block.derivative_block(1, 3 + m) =
          (pypt(1, m) * y[2] - y[1] * pypt(2, m)) / y[2] / y[2];
      }

      jacobian_matrix.camera_derivatives().push_back(camera_block);
      cameras_triplets.push_back(
        TripletType(j, k, jacobian_matrix.camera_derivatives().size()));

      typename JacobianMatrix::PointDerivativeBlock point_block;
      point_block.point_id = k;
      point_block.feature_id = i;
      for (Index m = 0; m < 3; m++)
      {
        point_block.derivative_block(0, m) =
          (pypp(0, m) * y[2] - y[0] * pypp(2, m)) / y[2] / y[2];
        point_block.derivative_block(1, m) =
          (pypp(1, m) * y[2] - y[1] * pypp(2, m)) / y[2] / y[2];
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

};

}//ba
}//sfm
}//hs

#endif
