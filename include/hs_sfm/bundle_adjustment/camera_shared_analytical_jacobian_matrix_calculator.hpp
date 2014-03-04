#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_ANALYTICAL_JACOBIAN_MATRIX_CALCULATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_CAMERA_SHARED_ANALYTICAL_JACOBIAN_MATRIX_CALCULATOR_HPP_

#include "hs_sfm/bundle_adjustment/camera_shared_jacobian_matrix.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _VectorFunction>
class CameraSharedAnalyticalJacobianMatrixCalculator;

template <typename _Scalar>
class CameraSharedAnalyticalJacobianMatrixCalculator<
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
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

  typedef EIGEN_MATRIX(Scalar,
                       VectorFunction::params_per_key_,
                       VectorFunction::params_per_point_)
          PointDerivative;
  typedef EIGEN_MATRIX(Scalar,
                       VectorFunction::params_per_key_,
                       VectorFunction::extrinsic_params_per_image_)
          ExtrinsicDerivative;
  typedef EIGEN_MATRIX(Scalar, VectorFunction::params_per_key_, 3)
          RadialDerivative;
  typedef EIGEN_MATRIX(Scalar, VectorFunction::params_per_key_, 2)
          DecenteringDerivative;
  typedef EIGEN_MATRIX(Scalar, VectorFunction::params_per_key_, 5)
          IntrinsicDerivative;


  typedef EIGEN_ROW_VECTOR(Scalar, VectorFunction::params_per_point_)
          RadiusPointDerivative;
  typedef EIGEN_ROW_VECTOR(Scalar,
                           VectorFunction::extrinsic_params_per_image_)
          RadiusExtrinsicDerivative;

public:
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
    Index x_point_begin = 0;
    Index x_image_begin =
      vector_function.is_fix_points() ? 0 : point_params_size;
    Index x_camera_begin =
      (vector_function.is_fix_points() ? 0 : point_params_size) +
      (vector_function.is_fix_images() ? 0 : extrinsic_params_size);
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

      Vector3 point =
        vector_function.is_fix_points() ?
        vector_function.fix_points()[point_id] :
        x.segment(point_id * VectorFunction::params_per_point_,
                  VectorFunction::params_per_point_);
      Vector3 rotation, translation;
      if (vector_function.is_fix_images())
      {
        rotation = vector_function.fix_images()[image_id].segment(0, 3);
        translation = vector_function.fix_images()[image_id].segment(3, 3);
      }
      else
      {
        rotation =
          x.segment(x_image_begin +
                   image_id * VectorFunction::extrinsic_params_per_image_, 3);
        translation =
          x.segment(x_image_begin +
                   image_id * VectorFunction::extrinsic_params_per_image_ + 3,
                   3);
      }

      VectorX intrinsic_params =
        vector_function.is_fix_cameras() ?
        vector_function.fix_cameras()[camera_id] :
        x.segment(x_camera_begin +
                  camera_id * intrinsic_params_size_per_camera,
                  intrinsic_params_size_per_camera);

      Vector2 normalized_key;
      VectorFunction::WorldPointToNormalizedKey(point,
                                                rotation,
                                                translation,
                                                normalized_key);

      PointDerivative pnpp;
      ExtrinsicDerivative pnpe;
      WorldPointToNormalizedPointDerivation(point, rotation, translation,
                                            pnpe, pnpp);

      Index intrinsic_offset = 0;
      RadialDerivative prdpr;
      PointDerivative prdpp;
      ExtrinsicDerivative prdpe;
      Vector2 rd;
      if (vector_function.intrinsic_computations_mask()[
            COMPUTE_RADIAL_DISTORTION])
      {
        Scalar k1 = intrinsic_params[intrinsic_offset + 0];
        Scalar k2 = intrinsic_params[intrinsic_offset + 1];
        Scalar k3 = intrinsic_params[intrinsic_offset + 2];
        VectorFunction::GetRadialDistortDelta(k1, k2, k3,
                                              normalized_key, rd);
        RadialDistortDerivation(normalized_key,
                                pnpp,
                                pnpe,
                                k1, k2, k3,
                                prdpp, prdpe, prdpr);
        intrinsic_offset += 3;
      }
      DecenteringDerivative pddpd;
      PointDerivative pddpp;
      ExtrinsicDerivative pddpe;
      Vector2 dd;
      if (vector_function.intrinsic_computations_mask()[
            COMPUTE_DECENTERING_DISTORTION])
      {
        Scalar d1 = intrinsic_params[intrinsic_offset + 0];
        Scalar d2 = intrinsic_params[intrinsic_offset + 1];
        VectorFunction::GetDecenteringDistortDelta(d1, d2,
                                                   normalized_key, dd);
        DecenteringDistortDerivation(normalized_key,
                                     pnpp,
                                     pnpe,
                                     d1, d2,
                                     pddpp, pddpe, pddpd);
        intrinsic_offset += 2;
      }

      if (vector_function.intrinsic_computations_mask()[
            COMPUTE_RADIAL_DISTORTION])
      {
        normalized_key += rd;
        pnpp += prdpp;
        pnpe += prdpe;
      }
      if (vector_function.intrinsic_computations_mask()[
            COMPUTE_DECENTERING_DISTORTION])
      {
        normalized_key += dd;
        pnpp += pddpp;
        pnpe += pddpe;
      }

      if (vector_function.intrinsic_computations_mask()[
            COMPUTE_INTRINSIC_PARAMS])
      {
        Scalar focal_length = intrinsic_params[intrinsic_offset + 0];
        Scalar skew = intrinsic_params[intrinsic_offset + 1];
        Scalar principal_x = intrinsic_params[intrinsic_offset + 2];
        Scalar principal_y = intrinsic_params[intrinsic_offset + 3];
        Scalar pixel_ratio = intrinsic_params[intrinsic_offset + 4];

        if (!vector_function.is_fix_points())
        {
          typename JacobianMatrix::PointDerivativeBlock point_block;
          point_block.point_id = point_id;
          point_block.key_id = i;
          IntrinsicPointDerivation(pnpp, focal_length, skew, pixel_ratio,
                                   point_block.derivative_block);
          jacobian_matrix.point_derivatives().push_back(point_block);
        }

        if (!vector_function.is_fix_images())
        {
          typename JacobianMatrix::ImageDerivativeBlock image_block;
          image_block.image_id = image_id;
          image_block.key_id = i;
          IntrinsicExtrinsicDerivation(pnpe, focal_length, skew, pixel_ratio,
                                       image_block.derivative_block);
          jacobian_matrix.image_derivatives().push_back(image_block);
        }

        if (!vector_function.is_fix_cameras())
        {
          typename JacobianMatrix::CameraDerivativeBlock camera_block;
          camera_block.derivative_block.resize(
            VectorFunction::params_per_key_, intrinsic_params_size_per_camera);
          camera_block.camera_id = camera_id;
          camera_block.key_id = i;
          intrinsic_offset = 0;
          if (vector_function.intrinsic_computations_mask()[
                COMPUTE_RADIAL_DISTORTION])
          {
            RadialDerivative pipr;
            IntrinsicRadialDerivation(prdpr, focal_length, skew, pixel_ratio,
                                      pipr);
            camera_block.derivative_block.block(0, intrinsic_offset,
                                                VectorFunction::params_per_key_,
                                                3) = pipr;
            intrinsic_offset += 3;
          }
          if (vector_function.intrinsic_computations_mask()[
                COMPUTE_DECENTERING_DISTORTION])
          {
            DecenteringDerivative pipd;
            IntrinsicDecenteringDerivation(
              pddpd, focal_length, skew, pixel_ratio, pipd);
            camera_block.derivative_block.block(0, intrinsic_offset,
                                                VectorFunction::params_per_key_,
                                                2) = pipd;
            intrinsic_offset += 2;
          }
          IntrinsicDerivative pipi;
          IntrinsicIntrinsicDerivation(normalized_key,
                                       focal_length,
                                       skew,
                                       principal_x,
                                       principal_y,
                                       pixel_ratio,
                                       pipi);
          camera_block.derivative_block.block(0, intrinsic_offset,
                                              VectorFunction::params_per_key_,
                                              5) = pipi;
          jacobian_matrix.camera_derivatives().push_back(camera_block);
        }
      }
      else
      {
        if (!vector_function.is_fix_points())
        {
          typename JacobianMatrix::PointDerivativeBlock point_block;
          point_block.point_id = point_id;
          point_block.key_id = i;
          point_block.derivative_block = pnpp;
          jacobian_matrix.point_derivatives().push_back(point_block);
        }

        if (!vector_function.is_fix_images())
        {
          typename JacobianMatrix::ImageDerivativeBlock image_block;
          image_block.image_id = image_id;
          image_block.key_id = i;
          image_block.derivative_block = pnpe;
          jacobian_matrix.image_derivatives().push_back(image_block);
        }

        if (!vector_function.is_fix_cameras())
        {
          typename JacobianMatrix::CameraDerivativeBlock camera_block;
          camera_block.derivative_block.resize(
            VectorFunction::params_per_key_, intrinsic_params_size_per_camera);
          camera_block.camera_id = camera_id;
          camera_block.key_id = i;
          intrinsic_offset = 0;
          if (vector_function.intrinsic_computations_mask()[
                COMPUTE_RADIAL_DISTORTION])
          {
            camera_block.derivative_block.block(0, intrinsic_offset,
                                                VectorFunction::params_per_key_,
                                                3) = prdpr;
            intrinsic_offset += 3;
          }
          if (vector_function.intrinsic_computations_mask()[
                COMPUTE_DECENTERING_DISTORTION])
          {
            camera_block.derivative_block.block(0, intrinsic_offset,
                                                VectorFunction::params_per_key_,
                                                2) = pddpd;
            intrinsic_offset += 2;
          }
          jacobian_matrix.camera_derivatives().push_back(camera_block);
        }
      }

    }// for (Index i = 0; i < number_of_keys; i++)

    if (!vector_function.is_fix_points())
      jacobian_matrix.SetPointConstraints(vector_function.point_constraints());
    if (!vector_function.is_fix_images())
      jacobian_matrix.SetImageConstraints(vector_function.image_constraints());
    if (!vector_function.is_fix_cameras())
      jacobian_matrix.SetCameraConstraints(
        vector_function.camera_constraints());

    return 0;
  }

  static Err WorldPointToNormalizedPointDerivation(
    const Vector3& p,
    const Vector3& r,
    const Vector3& t,
    ExtrinsicDerivative& extrinsic_derivative,
    PointDerivative& point_derivative)
  {
    Scalar theta = r.norm();
    Scalar theta2 = theta * theta;
    Scalar theta4 = theta2 * theta2;
    Scalar d = r.dot(p);

    //sin(theta)
    Scalar st = std::sin(theta);
    //cos(theta)
    Scalar ct = std::cos(theta);

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

    for (Index m = 0; m < 3; m++)
    {
      //for r
      extrinsic_derivative(0, m) =
        (pypr(0, m) * y[2] - y[0] * pypr(2, m)) / y[2] / y[2];
      extrinsic_derivative(1, m) =
        (pypr(1, m) * y[2] - y[1] * pypr(2, m)) / y[2] / y[2];

      //for t
      extrinsic_derivative(0, 3 + m) =
        (pypt(0, m) * y[2] - y[0] * pypt(2, m)) / y[2] / y[2];
      extrinsic_derivative(1, 3 + m) =
        (pypt(1, m) * y[2] - y[1] * pypt(2, m)) / y[2] / y[2];
    }

    for (Index m = 0; m < 3; m++)
    {
      point_derivative(0, m) =
        (pypp(0, m) * y[2] - y[0] * pypp(2, m)) / y[2] / y[2];
      point_derivative(1, m) =
        (pypp(1, m) * y[2] - y[1] * pypp(2, m)) / y[2] / y[2];
    }

    return 0;
  }

  static Err RadialDistortDerivation(const Vector2& n,
                                     const PointDerivative& pnpp,
                                     const ExtrinsicDerivative& pnpe,
                                     Scalar k1, Scalar k2, Scalar k3,
                                     PointDerivative& prdpp,
                                     ExtrinsicDerivative& prdpe,
                                     RadialDerivative& radial_derivative)
  {
    Scalar r2 = n.squaredNorm();
    Scalar r4 = r2 * r2;
    Scalar r6 = r2 * r4;

    RadiusPointDerivative pr2pp, pr4pp, pr6pp;
    RadiusExtrinsicDerivative pr2pe, pr4pe, pr6pe;

    pr2pp = Scalar(2) * n[0] * pnpp.row(0) + Scalar(2) * n[1] * pnpp.row(1);
    pr4pp = Scalar(2) * r2 * pr2pp;
    pr6pp = Scalar(3) * r4 * pr2pp;

    pr2pe = Scalar(2) * n[0] * pnpe.row(0) + Scalar(2) * n[1] * pnpe.row(1);
    pr4pe = Scalar(2) * r2 * pr2pe;
    pr6pe = Scalar(3) * r4 * pr2pe;

    Scalar coeff = k1 * r2 + k2 * r4 + k3 * r6;

    RadiusPointDerivative coeff_point_derivative =
      k1 * pr2pp + k2 * pr4pp + k3 * pr6pp;
    prdpp.row(0) = coeff_point_derivative * n[0] + coeff * pnpp.row(0);
    prdpp.row(1) = coeff_point_derivative * n[1] + coeff * pnpp.row(1);

    RadiusExtrinsicDerivative coeff_extrinsic_derivative =
      k1 * pr2pe + k2 * pr4pe + k3 * pr6pe;
    prdpe.row(0) = coeff_extrinsic_derivative * n[0] + coeff * pnpe.row(0);
    prdpe.row(1) = coeff_extrinsic_derivative * n[1] + coeff * pnpe.row(1);

    radial_derivative << r2 * n[0], r4 * n[0], r6 * n[0],
                         r2 * n[1], r4 * n[1], r6 * n[1];

    return 0;
  }

  static Err DecenteringDistortDerivation(
    const Vector2& n,
    const PointDerivative& pnpp,
    const ExtrinsicDerivative& pnpe,
    Scalar d1, Scalar d2,
    PointDerivative& pddpp,
    ExtrinsicDerivative& pddpe,
    DecenteringDerivative& decentering_derivative)
  {
    Scalar r2 = n.squaredNorm();

    RadiusPointDerivative pr2pp =
      Scalar(2) * n[0] * pnpp.row(0) + Scalar(2) * n[1] * pnpp.row(1);
    RadiusExtrinsicDerivative pr2pe =
      Scalar(2) * n[0] * pnpe.row(0) + Scalar(2) * n[1] * pnpe.row(1);

    pddpp.row(0) = Scalar(2) * d1 * (pnpp.row(0) * n[1] + pnpp.row(1) * n[0]) +
                   d2 * (pr2pp + Scalar(4) * n[0] * pnpp.row(0));
    pddpp.row(1) = Scalar(2) * d2 * (pnpp.row(0) * n[1] + pnpp.row(1) * n[0]) +
                   d1 * (pr2pp + Scalar(4) * n[1] * pnpp.row(1));

    pddpe.row(0) = Scalar(2) * d1 * (pnpe.row(0) * n[1] + pnpe.row(1) * n[0]) +
                   d2 * (pr2pe + Scalar(4) * n[0] * pnpe.row(0));
    pddpe.row(1) = Scalar(2) * d2 * (pnpe.row(0) * n[1] + pnpe.row(1) * n[0]) +
                   d1 * (pr2pe + Scalar(4) * n[1] * pnpe.row(1));

    decentering_derivative
      << Scalar(2) * n[0] * n[1], r2 + Scalar(2) * n[0] * n[0],
         r2 + Scalar(2) * n[1] * n[1], Scalar(2) * n[0] * n[1];

    return 0;
  }

  static void IntrinsicPointDerivation(const PointDerivative& pnpp,
                                      Scalar focal_length,
                                      Scalar skew,
                                      Scalar pixel_ratio,
                                      PointDerivative& pipp)
  {
    pipp.row(0) = focal_length * pnpp.row(0) + skew * pnpp.row(1);
    pipp.row(1) = focal_length * pixel_ratio * pnpp.row(1);
  }

  static void IntrinsicExtrinsicDerivation(const ExtrinsicDerivative& pnpe,
                                          Scalar focal_length,
                                          Scalar skew,
                                          Scalar pixel_ratio,
                                          ExtrinsicDerivative& pipe)
  {

    pipe.row(0) = focal_length * pnpe.row(0) + skew * pnpe.row(1);
    pipe.row(1) = focal_length * pixel_ratio * pnpe.row(1);
  }

  static void IntrinsicRadialDerivation(const RadialDerivative& prdpr,
                                       Scalar focal_length,
                                       Scalar skew,
                                       Scalar pixel_ratio,
                                       RadialDerivative& pipr)
  {
    pipr.row(0) = focal_length * prdpr.row(0) + skew * prdpr.row(1);
    pipr.row(1) = focal_length * pixel_ratio * prdpr.row(1);
  }

  static void IntrinsicDecenteringDerivation(const DecenteringDerivative& pddpd,
                                            Scalar focal_length,
                                            Scalar skew,
                                            Scalar pixel_ratio,
                                            DecenteringDerivative& pipd)
  {
    pipd.row(0) = focal_length * pddpd.row(0) + skew * pddpd.row(1);
    pipd.row(1) = focal_length * pixel_ratio * pddpd.row(1);
  }

  static void IntrinsicIntrinsicDerivation(const Vector2& n,
                                          Scalar focal_length,
                                          Scalar skew,
                                          Scalar principal_x,
                                          Scalar principal_y,
                                          Scalar pixel_ratio,
                                          IntrinsicDerivative& pipi)
  {
    pipi <<
      n[0], n[1], Scalar(1), Scalar(0), Scalar(0),
      pixel_ratio * n[1], Scalar(0), Scalar(0), Scalar(1), focal_length * n[1];
  }

};

}
}
}

#endif
