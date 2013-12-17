#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_GENERAL_SYNTHETIC_DATA_GENERATOR_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_GENERAL_SYNTHETIC_DATA_GENERATOR_HPP_

#include "hs_sfm/synthetic/scene_generator.hpp"
#include "hs_sfm/synthetic/keyset_generator.hpp"
#include "hs_sfm/bundle_adjustment/general_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar, typename _ImageDimension>
class GeneralSyntheticDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;

  typedef hs::sfm::synthetic::SceneGenerator<Scalar, ImageDimension>
          SceneGenerator;
  typedef typename SceneGenerator::IntrinsicParams IntrinsicParams;
  typedef typename SceneGenerator::IntrinsicParamsContainer
                   IntrinsicParamsContainer;
  typedef typename SceneGenerator::ExtrinsicParams ExtrinsicParams;
  typedef typename SceneGenerator::ExtrinsicParamsContainer
                   ExtrinsicParamsContainer;
  typedef typename SceneGenerator::Image Image;
  typedef typename SceneGenerator::ImageContainer ImageContainer;
  typedef typename SceneGenerator::Point3D Point3D;
  typedef typename SceneGenerator::Point3DContainer Point3DContainer;

  typedef hs::sfm::synthetic::KeysetGenerator<Scalar, ImageDimension>
          KeysetGenerator;
  typedef typename KeysetGenerator::Keyset Keyset;
  typedef typename KeysetGenerator::KeysetContainer KeysetContainer;

  typedef GeneralVectorFunction<Scalar> VectorFunction;
  typedef typename VectorFunction::Index Index;
  typedef typename VectorFunction::XVector XVector;
  typedef typename VectorFunction::YVector YVector;
  typedef typename VectorFunction::FeatureMap FeatureMap;
  typedef typename VectorFunction::FeatureMapContainer FeatureMapContainer;
  typedef typename VectorFunction::ImageCameraMap ImageCameraMap;
  typedef typename VectorFunction::IntrinsicComputationsMask
                   IntrinsicComputationsMask;
  typedef typename VectorFunction::IntrinsicConstraintsMask
                   IntrinsicConstraintsMask;
  typedef typename VectorFunction::StructureConstraintsMask
                   StructureConstraintsMask;

public:
  GeneralSyntheticDataGenerator() {}

private:
  SceneGenerator scene_generator;
  KeysetGenerator keyset_generator;
  size_t number_of_constrained_points_;
  size_t number_of_cameras_;
  IntrinsicComputationsMask intrinsic_computation_mask_;
  IntrinsicConstraintsMask intrinsic_constraints_mask_;
  StructureConstraintsMask structure_constraints_mask_;

};

}
}
}

#endif
