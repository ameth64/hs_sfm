#include <iostream>
#include <iomanip>

#include "gtest/gtest.h"

#include "hs_math/linear_algebra/latraits/vector_eigen.hpp"
#include "hs_math/linear_algebra/lafunc/arithmetic_eigen.hpp"
#include "hs_math/fdjac/forward_finite_difference_jacobian_matrix_calculator.hpp"

#include "hs_sfm/bundle_adjustment/camera_shared_analytical_jacobian_matrix_calculator.hpp"

namespace
{

template <typename _Scalar>
class TestCameraSharedAnalyticalDerivative
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef hs::sfm::ba::CameraSharedVectorFunction<Scalar>
          CameraSharedVectorFunction;
  typedef hs::sfm::ba::CameraSharedAnalyticalJacobianMatrixCalculator<
            CameraSharedVectorFunction>
          CameraSharedJacobianCalculator;
  typedef typename CameraSharedVectorFunction::Index Index;

  typedef EIGEN_MATRIX(Scalar,
                       CameraSharedVectorFunction::params_per_key_,
                       CameraSharedVectorFunction::extrinsic_params_per_image_)
          ExtrinsicDerivative;
  typedef EIGEN_MATRIX(Scalar,
                       CameraSharedVectorFunction::params_per_key_,
                       CameraSharedVectorFunction::params_per_point_)
          PointDerivative;
  typedef EIGEN_MATRIX(Scalar, CameraSharedVectorFunction::params_per_key_, 3)
          RadialDerivative;
  typedef EIGEN_MATRIX(Scalar, CameraSharedVectorFunction::params_per_key_, 2)
          DecenteringDerivative;
  typedef EIGEN_MATRIX(Scalar, CameraSharedVectorFunction::params_per_key_, 5)
          IntrinsicDerivative;

  typedef EIGEN_VECTOR(Scalar, 3) Vector3;
  typedef EIGEN_VECTOR(Scalar, 2) Vector2;
  typedef EIGEN_VECTOR(Scalar, Eigen::Dynamic) VectorX;
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

  class NormalizedKeyVectorFunction
  {
  public:
    typedef int Err;
    typedef EIGEN_VECTOR(
              Scalar,
              CameraSharedVectorFunction::params_per_point_ +
              CameraSharedVectorFunction::extrinsic_params_per_image_)
            XVector;
    typedef EIGEN_VECTOR(Scalar,
                         CameraSharedVectorFunction::params_per_key_)
            YVector;

    Err operator() (const XVector& x, YVector& y) const
    {
      Vector3 point =
        x.segment(0, CameraSharedVectorFunction::params_per_point_);
      Vector3 rotation =
        x.segment(CameraSharedVectorFunction::params_per_point_, 3);
      Vector3 translation =
        x.segment(CameraSharedVectorFunction::params_per_point_ + 3, 3);
      CameraSharedVectorFunction::WorldPointToNormalizedKey(
        point,
        rotation,
        translation,
        y);
      return 0;
    }
  };
  typedef hs::math::fdjac::ForwardFiniteDifferenceJacobianMatrixCalculator<
            NormalizedKeyVectorFunction>
          NormalizedKeyFFDJacobianCalculator;

  class NormalizedKeyAnalyticalJacobianCalculator
  {
  public:
    typedef NormalizedKeyVectorFunction VectorFunction;
    typedef typename VectorFunction::XVector XVector;
    typedef EIGEN_MATRIX(
              Scalar,
              CameraSharedVectorFunction::params_per_key_,
              CameraSharedVectorFunction::extrinsic_params_per_image_ +
              CameraSharedVectorFunction::params_per_point_)
            JacobianMatrix;

    Err operator() (const VectorFunction& vector_function,
                    const XVector& x,
                    JacobianMatrix& jacobian_matrix) const
    {
      Vector3 point =
        x.segment(0, CameraSharedVectorFunction::params_per_point_);
      Vector3 rotation =
        x.segment(CameraSharedVectorFunction::params_per_point_, 3);
      Vector3 translation =
        x.segment(CameraSharedVectorFunction::params_per_point_ + 3, 3);
      ExtrinsicDerivative extrinsic_derivative;
      PointDerivative point_derivative;
      CameraSharedJacobianCalculator::WorldPointToNormalizedPointDerivation(
        point, rotation, translation,
        extrinsic_derivative, point_derivative);
      jacobian_matrix.block(
        0, 0,
        CameraSharedVectorFunction::params_per_key_,
        CameraSharedVectorFunction::params_per_point_) = point_derivative;
      jacobian_matrix.block(
        0, CameraSharedVectorFunction::params_per_point_,
        CameraSharedVectorFunction::params_per_key_,
        CameraSharedVectorFunction::extrinsic_params_per_image_) =
        extrinsic_derivative;
      return 0;
    }
  };

  class RadialVectorFunction
  {
  public:
    typedef int Err;

    typedef EIGEN_VECTOR(
              Scalar,
              CameraSharedVectorFunction::params_per_point_ +
              CameraSharedVectorFunction::extrinsic_params_per_image_ + 3)
            XVector;
    typedef EIGEN_VECTOR(Scalar,
                         CameraSharedVectorFunction::params_per_key_)
            YVector;

    Err operator() (const XVector& x, YVector& y) const
    {
      Vector3 point =
        x.segment(0, CameraSharedVectorFunction::params_per_point_);
      Vector3 rotation =
        x.segment(CameraSharedVectorFunction::params_per_point_, 3);
      Vector3 translation =
        x.segment(CameraSharedVectorFunction::params_per_point_ + 3, 3);
      Vector2 normalized_key;
      CameraSharedVectorFunction::WorldPointToNormalizedKey(
        point,
        rotation,
        translation,
        normalized_key);

      Vector2 delta;
      CameraSharedVectorFunction::GetRadialDistortDelta(
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 0],
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 1],
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 2],
        normalized_key,
        delta);

      y = normalized_key + delta;

      return 0;
    }
  };
  typedef hs::math::fdjac::ForwardFiniteDifferenceJacobianMatrixCalculator<
            RadialVectorFunction>
          RadialFFDJacobianCalculator;

  class RadialAnalyticalJacobianCalculator
  {
  public:
    typedef RadialVectorFunction VectorFunction;
    typedef typename VectorFunction::XVector XVector;
    typedef EIGEN_MATRIX(
              Scalar,
              CameraSharedVectorFunction::params_per_key_,
              CameraSharedVectorFunction::params_per_point_ +
              CameraSharedVectorFunction::extrinsic_params_per_image_ + 3)
            JacobianMatrix;

    Err operator() (const VectorFunction& vector_function,
                    const XVector& x,
                    JacobianMatrix& jacobian_matrix) const
    {
      typename NormalizedKeyVectorFunction::XVector normal_x;
      normal_x = x.segment(0, normal_x.rows());
      NormalizedKeyVectorFunction normal_vector_function;
      typename NormalizedKeyVectorFunction::YVector normal_y;
      if (normal_vector_function(normal_x, normal_y) != 0)
      {
        return -1;
      }
      NormalizedKeyAnalyticalJacobianCalculator normal_calculator;
      typename NormalizedKeyAnalyticalJacobianCalculator::JacobianMatrix
        normal_jacobian_matrix;
      if (normal_calculator(normal_vector_function,
                            normal_x,
                            normal_jacobian_matrix) != 0)
      {
        return -1;
      }

      Scalar k1 = x[normal_x.rows() + 0];
      Scalar k2 = x[normal_x.rows() + 1];
      Scalar k3 = x[normal_x.rows() + 2];

      PointDerivative pnpp =
        normal_jacobian_matrix.block(
          0, 0,
          CameraSharedVectorFunction::params_per_key_,
          CameraSharedVectorFunction::params_per_point_);
      ExtrinsicDerivative pnpe =
        normal_jacobian_matrix.block(
          0, CameraSharedVectorFunction::params_per_point_,
          CameraSharedVectorFunction::params_per_key_,
          CameraSharedVectorFunction::extrinsic_params_per_image_);

      PointDerivative prdpp;
      ExtrinsicDerivative prdpe;
      RadialDerivative radial_derivative;

      CameraSharedJacobianCalculator::RadialDistortDerivation(
        normal_y, pnpp, pnpe, k1, k2, k3, prdpp, prdpe, radial_derivative);

      jacobian_matrix.block(
        0, 0,
        CameraSharedVectorFunction::params_per_key_,
        CameraSharedVectorFunction::params_per_point_) = pnpp + prdpp;
      jacobian_matrix.block(
        0, CameraSharedVectorFunction::params_per_point_,
        CameraSharedVectorFunction::params_per_key_,
        CameraSharedVectorFunction::extrinsic_params_per_image_) = pnpe + prdpe;
      jacobian_matrix.block(
        0,
        CameraSharedVectorFunction::params_per_point_ +
        CameraSharedVectorFunction::extrinsic_params_per_image_,
        CameraSharedVectorFunction::params_per_key_, 3) = radial_derivative;

      return 0;
    }
  };

  class DecenteringVectorFunction
  {
  public:
    typedef int Err;
    typedef EIGEN_VECTOR(
              Scalar,
              CameraSharedVectorFunction::params_per_point_ +
              CameraSharedVectorFunction::extrinsic_params_per_image_ + 2)
            XVector;
    typedef EIGEN_VECTOR(Scalar,
                         CameraSharedVectorFunction::params_per_key_)
            YVector;

    Err operator() (const XVector& x, YVector& y) const
    {
      Vector3 point =
        x.segment(0, CameraSharedVectorFunction::params_per_point_);
      Vector3 rotation =
        x.segment(CameraSharedVectorFunction::params_per_point_, 3);
      Vector3 translation =
        x.segment(CameraSharedVectorFunction::params_per_point_ + 3, 3);
      Vector2 normalized_key;
      CameraSharedVectorFunction::WorldPointToNormalizedKey(
        point,
        rotation,
        translation,
        normalized_key);

      Vector2 delta;
      CameraSharedVectorFunction::GetDecenteringDistortDelta(
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 0],
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 1],
        normalized_key,
        delta);

      y = normalized_key + delta;

      return 0;
    }
  };
  typedef hs::math::fdjac::ForwardFiniteDifferenceJacobianMatrixCalculator<
            DecenteringVectorFunction>
          DecenteringFFDJacobianCalculator;

  class DecenteringAnalyticalJacobianCalculator
  {
  public:
    typedef DecenteringVectorFunction VectorFunction;
    typedef typename VectorFunction::XVector XVector;
    typedef EIGEN_MATRIX(
              Scalar,
              CameraSharedVectorFunction::params_per_key_,
              CameraSharedVectorFunction::params_per_point_ +
              CameraSharedVectorFunction::extrinsic_params_per_image_ + 2)
            JacobianMatrix;

    Err operator() (const VectorFunction& vector_function,
                    const XVector& x,
                    JacobianMatrix& jacobian_matrix) const
    {
      typename NormalizedKeyVectorFunction::XVector normal_x;
      normal_x = x.segment(0, normal_x.rows());
      NormalizedKeyVectorFunction normal_vector_function;
      typename NormalizedKeyVectorFunction::YVector normal_y;
      if (normal_vector_function(normal_x, normal_y) != 0)
      {
        return -1;
      }
      NormalizedKeyAnalyticalJacobianCalculator normal_calculator;
      typename NormalizedKeyAnalyticalJacobianCalculator::JacobianMatrix
        normal_jacobian_matrix;
      if (normal_calculator(normal_vector_function,
                            normal_x,
                            normal_jacobian_matrix) != 0)
      {
        return -1;
      }

      Scalar d1 = x[normal_x.rows() + 0];
      Scalar d2 = x[normal_x.rows() + 1];

      PointDerivative pnpp =
        normal_jacobian_matrix.block(
          0, 0,
          CameraSharedVectorFunction::params_per_key_,
          CameraSharedVectorFunction::params_per_point_);
      ExtrinsicDerivative pnpe =
        normal_jacobian_matrix.block(
          0, CameraSharedVectorFunction::params_per_point_,
          CameraSharedVectorFunction::params_per_key_,
          CameraSharedVectorFunction::extrinsic_params_per_image_);

      PointDerivative pddpp;
      ExtrinsicDerivative pddpe;
      DecenteringDerivative decentering_derivative;

      CameraSharedJacobianCalculator::DecenteringDistortDerivation(
        normal_y, pnpp, pnpe, d1, d2, pddpp, pddpe, decentering_derivative);

      jacobian_matrix.block(
        0, 0,
        CameraSharedVectorFunction::params_per_key_,
        CameraSharedVectorFunction::params_per_point_) = pnpp + pddpp;
      jacobian_matrix.block(
        0, CameraSharedVectorFunction::params_per_point_,
        CameraSharedVectorFunction::params_per_key_,
        CameraSharedVectorFunction::extrinsic_params_per_image_) = pnpe + pddpe;
      jacobian_matrix.block(
        0,
        CameraSharedVectorFunction::params_per_point_ +
        CameraSharedVectorFunction::extrinsic_params_per_image_,
        CameraSharedVectorFunction::params_per_key_, 2) =
        decentering_derivative;

      return 0;
    }
  };

  class IntrinsicVectorFunction
  {
  public:
    typedef int Err;

    typedef EIGEN_VECTOR(
              Scalar,
              CameraSharedVectorFunction::params_per_point_ +
              CameraSharedVectorFunction::extrinsic_params_per_image_ + 10)
            XVector;
    typedef EIGEN_VECTOR(Scalar, CameraSharedVectorFunction::params_per_key_)
            YVector;

    Err operator() (const XVector& x, YVector& y) const
    {
      Vector3 point =
        x.segment(0, CameraSharedVectorFunction::params_per_point_);
      Vector3 rotation =
        x.segment(CameraSharedVectorFunction::params_per_point_, 3);
      Vector3 translation =
        x.segment(CameraSharedVectorFunction::params_per_point_ + 3, 3);
      Vector2 normalized_key;
      CameraSharedVectorFunction::WorldPointToNormalizedKey(
        point,
        rotation,
        translation,
        normalized_key);

      Vector2 radial_delta;
      CameraSharedVectorFunction::GetRadialDistortDelta(
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 0],
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 1],
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 2],
        normalized_key,
        radial_delta);

      Vector2 decentering_delta;
      CameraSharedVectorFunction::GetDecenteringDistortDelta(
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 3],
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 4],
        normalized_key,
        decentering_delta);

      normalized_key += radial_delta + decentering_delta;

      CameraSharedVectorFunction::NormalizedKeyToImageKey(
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 5],
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 6],
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 7],
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 8],
        x[CameraSharedVectorFunction::params_per_point_ +
          CameraSharedVectorFunction::extrinsic_params_per_image_ + 9],
        normalized_key, y);

      return 0;
    }
  };
  typedef hs::math::fdjac::ForwardFiniteDifferenceJacobianMatrixCalculator<
            IntrinsicVectorFunction>
          IntrinsicFFDJacobianCalculator;

  class IntrinsicAnalyticalJacobianCalculator
  {
  public:
    typedef IntrinsicVectorFunction VectorFunction;
    typedef typename VectorFunction::XVector XVector;
    typedef EIGEN_MATRIX(
              Scalar,
              CameraSharedVectorFunction::params_per_key_,
              CameraSharedVectorFunction::params_per_point_ +
              CameraSharedVectorFunction::extrinsic_params_per_image_ + 10)
            JacobianMatrix;

    Err operator() (const VectorFunction& vector_function,
                    const XVector& x,
                    JacobianMatrix& jacobian_matrix) const
    {
      typename NormalizedKeyVectorFunction::XVector normal_x;
      normal_x = x.segment(0, normal_x.rows());
      NormalizedKeyVectorFunction normal_vector_function;
      typename NormalizedKeyVectorFunction::YVector normal_y;
      if (normal_vector_function(normal_x, normal_y) != 0)
      {
        return -1;
      }
      NormalizedKeyAnalyticalJacobianCalculator normal_calculator;
      typename NormalizedKeyAnalyticalJacobianCalculator::JacobianMatrix
        normal_jacobian_matrix;
      if (normal_calculator(normal_vector_function,
                            normal_x,
                            normal_jacobian_matrix) != 0)
      {
        return -1;
      }

      PointDerivative pnpp =
        normal_jacobian_matrix.block(
          0, 0,
          CameraSharedVectorFunction::params_per_key_,
          CameraSharedVectorFunction::params_per_point_);
      ExtrinsicDerivative pnpe =
        normal_jacobian_matrix.block(
          0, CameraSharedVectorFunction::params_per_point_,
          CameraSharedVectorFunction::params_per_key_,
          CameraSharedVectorFunction::extrinsic_params_per_image_);

      Scalar k1 = x[normal_x.rows() + 0];
      Scalar k2 = x[normal_x.rows() + 1];
      Scalar k3 = x[normal_x.rows() + 2];
      Scalar d1 = x[normal_x.rows() + 3];
      Scalar d2 = x[normal_x.rows() + 4];
      Scalar focal_length = x[normal_x.rows() + 5];
      Scalar skew = x[normal_x.rows() + 6];
      Scalar principal_x = x[normal_x.rows() + 7];
      Scalar principal_y = x[normal_x.rows() + 8];
      Scalar pixel_ratio = x[normal_x.rows() + 9];

      Vector2 rd;
      CameraSharedVectorFunction::GetRadialDistortDelta(
        k1, k2, k3, normal_y, rd);

      Vector2 dd;
      CameraSharedVectorFunction::GetDecenteringDistortDelta(
        d1, d2, normal_y, dd);

      PointDerivative prdpp, pddpp;
      ExtrinsicDerivative prdpe, pddpe;
      RadialDerivative prdpr;
      DecenteringDerivative pddpd;

      if (CameraSharedJacobianCalculator::RadialDistortDerivation(
            normal_y, pnpp, pnpe, k1, k2, k3,
            prdpp, prdpe, prdpr) != 0)
      {
        return -1;
      }

      if (CameraSharedJacobianCalculator::DecenteringDistortDerivation(
            normal_y, pnpp, pnpe, d1, d2,
            pddpp, pddpe, pddpd) != 0)
      {
        return -1;
      }

      PointDerivative pipp;
      ExtrinsicDerivative pipe;
      RadialDerivative pipr;
      DecenteringDerivative pipd;
      IntrinsicDerivative pipi;

      pnpp += prdpp + pddpp;
      pnpe += prdpe + pddpe;
      CameraSharedJacobianCalculator::IntrinsicPointDerivation(
        pnpp, focal_length, skew, pixel_ratio, pipp);
      CameraSharedJacobianCalculator::IntrinsicExtrinsicDerivation(
        pnpe, focal_length, skew, pixel_ratio, pipe);
      CameraSharedJacobianCalculator::IntrinsicRadialDerivation(
        prdpr, focal_length, skew, pixel_ratio, pipr);
      CameraSharedJacobianCalculator::IntrinsicDecenteringDerivation(
        pddpd, focal_length, skew, pixel_ratio, pipd);
      normal_y += rd + dd;
      CameraSharedJacobianCalculator::IntrinsicIntrinsicDerivation(
        normal_y, focal_length, skew, principal_x, principal_y, pixel_ratio,
        pipi);

      jacobian_matrix.block(
        0, 0,
        CameraSharedVectorFunction::params_per_key_,
        CameraSharedVectorFunction::params_per_point_) = pipp;
      jacobian_matrix.block(
        0, CameraSharedVectorFunction::params_per_point_,
        CameraSharedVectorFunction::params_per_key_,
        CameraSharedVectorFunction::extrinsic_params_per_image_) = pipe;
      jacobian_matrix.block(
        0,
        CameraSharedVectorFunction::params_per_point_ +
        CameraSharedVectorFunction::extrinsic_params_per_image_,
        CameraSharedVectorFunction::params_per_key_, 3) = pipr;
      jacobian_matrix.block(
        0,
        CameraSharedVectorFunction::params_per_point_ +
        CameraSharedVectorFunction::extrinsic_params_per_image_ + 3,
        CameraSharedVectorFunction::params_per_key_, 2) = pipd;
      jacobian_matrix.block(
        0,
        CameraSharedVectorFunction::params_per_point_ +
        CameraSharedVectorFunction::extrinsic_params_per_image_ + 5,
        CameraSharedVectorFunction::params_per_key_, 5) = pipi;

      return 0;
    }
  };

public:
  Err TestNormalizedKey(const Vector3& point,
                        const Vector3& rotation,
                        const Vector3& translation) const
  {
    NormalizedKeyVectorFunction vector_function;
    typename NormalizedKeyVectorFunction::XVector x;
    x.segment(0, CameraSharedVectorFunction::params_per_point_) = point;
    x.segment(CameraSharedVectorFunction::params_per_point_, 3) = rotation;
    x.segment(CameraSharedVectorFunction::params_per_point_ + 3, 3) =
      translation;

    NormalizedKeyFFDJacobianCalculator ffd_calculator(Scalar(1e-6),
                                                      Scalar(1e-8));
    NormalizedKeyAnalyticalJacobianCalculator analytical_calculator;

    typename NormalizedKeyFFDJacobianCalculator::JacobianMatrix
      ffd_jacobian_matrix;
    typename NormalizedKeyAnalyticalJacobianCalculator::JacobianMatrix
      analytical_jacobian_matrix;

    if (ffd_calculator(vector_function, x, ffd_jacobian_matrix) != 0)
      return -1;
    if (analytical_calculator(vector_function, x,
                              analytical_jacobian_matrix) != 0)
      return -1;

    if (!ffd_jacobian_matrix.isApprox(analytical_jacobian_matrix, Scalar(1e-5)))
      return -1;

    return 0;
  }

  Err TestRadial(const Vector3& point,
                 const Vector3& rotation,
                 const Vector3& translation,
                 Scalar k1, Scalar k2, Scalar k3) const
  {
    RadialVectorFunction vector_function;
    typename RadialVectorFunction::XVector x;
    x.segment(0, CameraSharedVectorFunction::params_per_point_) = point;
    x.segment(CameraSharedVectorFunction::params_per_point_, 3) = rotation;
    x.segment(CameraSharedVectorFunction::params_per_point_ + 3, 3) =
      translation;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 0] = k1;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 1] = k2;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 2] = k3;

    RadialFFDJacobianCalculator ffd_calculator(Scalar(1e-6), Scalar(1e-8));
    RadialAnalyticalJacobianCalculator analytical_calculator;

    typename RadialFFDJacobianCalculator::JacobianMatrix
      ffd_jacobian_matrix;
    typename RadialAnalyticalJacobianCalculator::JacobianMatrix
      analytical_jacobian_matrix;

    if (ffd_calculator(vector_function, x, ffd_jacobian_matrix) != 0)
      return -1;
    if (analytical_calculator(vector_function, x,
                              analytical_jacobian_matrix) != 0)
      return -1;

    if (!ffd_jacobian_matrix.isApprox(analytical_jacobian_matrix, Scalar(1e-5)))
      return -1;

    return 0;
  }

  Err TestDecentering(const Vector3& point,
                      const Vector3& rotation,
                      const Vector3& translation,
                      Scalar d1, Scalar d2) const
  {
    DecenteringVectorFunction vector_function;
    typename DecenteringVectorFunction::XVector x;
    x.segment(0, CameraSharedVectorFunction::params_per_point_) = point;
    x.segment(CameraSharedVectorFunction::params_per_point_, 3) = rotation;
    x.segment(CameraSharedVectorFunction::params_per_point_ + 3, 3) =
      translation;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 0] = d1;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 1] = d2;

    DecenteringFFDJacobianCalculator ffd_calculator(Scalar(1e-6), Scalar(1e-8));
    DecenteringAnalyticalJacobianCalculator analytical_calculator;

    typename DecenteringFFDJacobianCalculator::JacobianMatrix
      ffd_jacobian_matrix;
    typename DecenteringAnalyticalJacobianCalculator::JacobianMatrix
      analytical_jacobian_matrix;

    if (ffd_calculator(vector_function, x, ffd_jacobian_matrix) != 0)
      return -1;
    if (analytical_calculator(vector_function, x,
                              analytical_jacobian_matrix) != 0)
      return -1;

    if (!ffd_jacobian_matrix.isApprox(analytical_jacobian_matrix, Scalar(1e-5)))
      return -1;

    return 0;
  }

  Err TestIntrinsic(const Vector3& point,
                    const Vector3& rotation,
                    const Vector3& translation,
                    Scalar k1, Scalar k2, Scalar k3,
                    Scalar d1, Scalar d2,
                    Scalar focal_length,
                    Scalar skew,
                    Scalar principal_x,
                    Scalar principal_y,
                    Scalar pixel_radio) const
  {
    IntrinsicVectorFunction vector_function;
    typename IntrinsicVectorFunction::XVector x;
    x.segment(0, CameraSharedVectorFunction::params_per_point_) = point;
    x.segment(CameraSharedVectorFunction::params_per_point_, 3) = rotation;
    x.segment(CameraSharedVectorFunction::params_per_point_ + 3, 3) =
      translation;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 0] = k1;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 1] = k2;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 2] = k3;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 3] = d1;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 4] = d2;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 5] =
      focal_length;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 6] =
      skew;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 7] =
      principal_x;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 8] =
      principal_y;
    x[CameraSharedVectorFunction::params_per_point_ +
      CameraSharedVectorFunction::extrinsic_params_per_image_ + 9] =
      pixel_radio;

    IntrinsicFFDJacobianCalculator ffd_calculator(Scalar(1e-6), Scalar(1e-8));
    IntrinsicAnalyticalJacobianCalculator analytical_calculator;

    typename IntrinsicFFDJacobianCalculator::JacobianMatrix
      ffd_jacobian_matrix;
    typename IntrinsicAnalyticalJacobianCalculator::JacobianMatrix
      analytical_jacobian_matrix;

    if (ffd_calculator(vector_function, x, ffd_jacobian_matrix) != 0)
      return -1;
    if (analytical_calculator(vector_function, x,
                              analytical_jacobian_matrix) != 0)
      return -1;

    std::cout.setf(std::ios::fixed);
    std::cout<<std::setprecision(10);
    std::cout<<"       ffd:"<<ffd_jacobian_matrix(0, 3)<<"\n";
    std::cout<<"analytical:"<<analytical_jacobian_matrix(0, 3)<<"\n";

    if (!ffd_jacobian_matrix.isApprox(analytical_jacobian_matrix, Scalar(1e-5)))
      return -1;

    return 0;
  }
};

TEST(TestCameraSharedAnalyticalDerivative, SimpleTest)
{
  typedef double Scalar;
  typedef TestCameraSharedAnalyticalDerivative<Scalar> Tester;

  Tester tester;

  Tester::Vector3 point;
  Tester::Vector3 rotation;
  Tester::Vector3 translation;

  point << 89.06332903711264,
           -254.7119654139254,
           27.36102981725708;
  rotation << 0.0461972342123952,
              -0.01398317345699249,
              -1.02689361003788;
  translation << 228.8962690432249,
                 77.85673263411233,
                 -476.5385356161123;

  Scalar k1 = -0.005049046200604883;
  Scalar k2 = 0.03129380429860988;
  Scalar k3 = -0.03079482096044222;

  Scalar d1 = -0.0005537654832018913;
  Scalar d2 = -0.0004987771776838148;

  Scalar focal_length = 4835.47665904517;
  Scalar skew = 0;
  Scalar principal_x = -42.4095312016;
  Scalar principal_y = -31.699212823;
  Scalar pixel_ratio = 1;

  ASSERT_EQ(0, tester.TestNormalizedKey(point, rotation, translation));
  ASSERT_EQ(0, tester.TestRadial(point, rotation, translation, k1, k2, k3));
  ASSERT_EQ(0, tester.TestDecentering(point, rotation, translation, d1, d2));
  ASSERT_EQ(0, tester.TestIntrinsic(point, rotation, translation,
                                    k1, k2, k3,
                                    d1, d2,
                                    focal_length,
                                    skew,
                                    principal_x,
                                    principal_y,
                                    pixel_ratio));
}

}
