#ifndef _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_VECTOR_FUNCTION_HPP_
#define _HS_SFM_BUNDLE_ADJUSTMENT_BA_GCP_CONSTRAINED_VECTOR_FUNCTION_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/bundle_adjustment/ba_naive_vector_function.hpp"

namespace hs
{
namespace sfm
{
namespace ba
{

template <typename _Scalar>
class BAGCPConstrainedVectorFunction : public BANaiveVectorFunction<_Scalar>
{
public:
  typedef _Scalar Scalar;
  typedef int Err;
  typedef BANaiveVectorFunction<Scalar> Base;
  typedef typename Base::Index Index;
  typedef typename Base::XVector XVector;
  typedef typename Base::YVector YVector;
  typedef typename Base::FeatureMap FeatureMap;
  typedef typename Base::FeatureMapContainer FeatureMapContainer;
  static const Index params_per_camera_ = Base::params_per_camera_;
  static const Index params_per_point_ = Base::params_per_point_;
  static const Index params_per_feature_ = Base::params_per_feature_;
  using Base::GetCameraParamsSize;
  using Base::GetXSize;
  using Base::GetYSize;
  using Base::feature_maps;
  using Base::number_of_features;
  using Base::PointProjectToFeature;

private:
  typedef typename Base::Vector2 Vector2;
  typedef typename Base::Vector3 Vector3;

public:
  BAGCPConstrainedVectorFunction() {}
  BAGCPConstrainedVectorFunction(Index number_of_cameras,
                                 Index number_of_points,
                                 Index number_of_features,
                                 const FeatureMapContainer& feature_maps,
                                 Index number_of_gcps)
    : Base(number_of_cameras, number_of_points, number_of_features,
           feature_maps),
      number_of_gcps_(number_of_gcps) {}

  inline Index number_of_gcps() const
  {
    return number_of_gcps_;
  }

  inline void set_number_of_gcps(Index number_of_gcps)
  {
    number_of_gcps_ = number_of_gcps;
  }

  //将像控点坐标作为测量值放在y向量的最后几个分量
  inline Index GetYSize() const
  {
    return Base::GetYSize() + number_of_gcps_ * Base::params_per_point_;
  }

  inline Index GetYFeatureSize() const
  {
    return Base::GetYSize();
  }

  Err operator() (const XVector& x, YVector& y) const
  {
    Index x_size = x.rows();
    if (x_size != GetXSize()) return -1;
    if (number_of_features() != Index(feature_maps().size())) return -1;
    Index y_size = GetYSize();
    y.resize(y_size);

    Index camera_params_size = GetCameraParamsSize();

    for (Index i = 0; i < number_of_features(); i++)
    {
      Index camera_id = feature_maps()[i].first;
      Index point_id = feature_maps()[i].second;

      Vector3 p = x.segment(camera_params_size + point_id * params_per_point_,
        params_per_point_);

      Vector3 r = x.segment(camera_id * params_per_camera_, 3);
      Vector3 t = x.segment(camera_id * params_per_camera_ + 3, 3);

      Vector2 feat;
      PointProjectToFeature(r, t, p, feat);

      y.segment(i * params_per_feature_, params_per_feature_) = feat;
    }

    Index y_gcp_start = y_size - number_of_gcps_ * Base::params_per_point_;
    Index x_gcp_start = x_size - number_of_gcps_ * Base::params_per_point_;

    for (Index i = 0; i < number_of_gcps_; i++)
    {
      y.segment(y_gcp_start + i * Base::params_per_point_,
                Base::params_per_point_) =
        x.segment(x_gcp_start + i * Base::params_per_point_,
                  Base::params_per_point_);
    }

    return 0;
  }

private:
  Index number_of_gcps_;
};

}
}
}

#endif
