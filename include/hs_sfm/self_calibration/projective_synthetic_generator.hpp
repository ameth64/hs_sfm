#ifndef _HS_SFM_SELF_CALIBRATION_PROJECTIVE_SYNTHETIC_GENERATOR_HPP_
#define _HS_SFM_SELF_CALIBRATION_PROJECTIVE_SYNTHETIC_GENERATOR_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"
#include "hs_sfm/synthetic/flight_generator.hpp"

namespace hs
{
namespace sfm
{
namespace sc
{

template <typename _Scalar, typename _ImageDimension>
class ProjectiveSyntheticGenrator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

  typedef hs::sfm::synthetic::FlightGenerator<Scalar, ImageDimension>
          FlightGenerator;
  typedef typename FlightGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename FlightGenerator::ExtrinsicParamsContainer
          ExtrinsicParamsContainer;
  typedef typename FlightGenerator::Image Image;
  typedef typename FlightGenerator::ImageContainer ImageContainer;
  typedef typename FlightGenerator::Point3D Point3D;
  typedef typename FlightGenerator::Point3DContainer Point3DContainer;

  typedef EIGEN_MATRIX(Scalar, 3, 4) PMatrix;
  typedef EIGEN_STD_VECTOR(PMatrix) PMatrixContainer;
};

}
}
}

#endif
