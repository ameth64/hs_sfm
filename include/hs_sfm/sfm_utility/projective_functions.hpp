#ifndef _HS_SFM_SFM_UTILITY_PROJECTIVE_FUNCTIONS_HPP_
#define _HS_SFM_SFM_UTILITY_PROJECTIVE_FUNCTIONS_HPP_

#include "hs_math/linear_algebra/eigen_macro.hpp"

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/radial_distortor.hpp"
#include "hs_sfm/sfm_utility/decentering_distortor.hpp"

namespace hs
{
namespace sfm
{

/**
 *  封装世界座标经过投影变换后映射到影像座标的函数。
 *
 *  世界座标到影像座标的投影变换可表示为：
 *  \f[
 *    \mathbf{x_i}=I(D(N(E(\mathbf{x_w}))))
 *  \f]
 *
 *  其中，\f$\mathbf{x_w}\f$是三维列向量，表示世界座标系的座标点。
 *
 *  \f$\mathbf{x_i}\f$是二维列向量，表示影像座标系的座标点，单位为像素，
 *  影像座标系原点为照片左上角像素的中心，x轴朝右，y轴朝下。
 *
 *  \f$E:R^3\mapsto R^2\f$
 *  表示由世界座标系的三维座标点到相机座标系的三维座标点的变换，
 *  该变换由[相机外参数](@ref CameraExtrinsicParams)确定。
 *
 *  \f$N:R^3\mapsto R^2\f$表示投影变换：
 *  \f[
 *    N(u, v, w) = \left[
 *    \begin{array}{l}
 *      u / w \\
 *      v / w
 *    \end{array}
 *    \right]
 *  \f]
 *
 *  \f$D:R^2\mapsto R^2\f$
 *  表示镜头畸变，
 *  包含[径向畸变](@ref RadialDistortor)和[偏心畸变](@ref DecenteringDistortor)
 *  两个部分，即
 *  \f[
 *    D(u, v) = \left[
 *      \begin{array}{l}
 *        u + (k_1r^2+k_2r^4+k_3r^6)u + 2d_1uv + d_2(r^2 + 2u^2) \\
 *        v + (k_1r^2+k_2r^4+k_3r^6)v + 2d_2uv + d_1(r^2 + 2v^2)
 *      \end{array}
 *    \right]
 *  \f]
 *  其中\f$r^2 = u^2 + v^2\f$，\f$k_1\f$、\f$k_2\f$、\f$k_3\f$是径向畸变参数，
 *  \f$d_1\f$、\f$d_2\f$是偏心畸变参数。
 *
 *  \f$I:R^2\mapsto R^2\f$
 *  表示相机内参数中的线性部分，即
 *  \f[
 *    I(u, v) = \left[
 *    \begin{array}{l}
 *      u * f + s * v + x \\
 *      f * r * v + y
 *    \end{array}
 *    \right]
 *  \f]
 *  其中\f$f\f$、\f$s\f$、\f$x\f$、\f$y\f$、\f$r\f$
 *  分别是[相机内参数](@ref CameraIntrinsicParams)中的
 *  焦距、偏斜、像主点x座标、像主点y座标和像素比率。
 *
 */
template <typename _Scalar>
class ProjectiveFunctions
{
public:
  typedef _Scalar Scalar;
  typedef int Err;

  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_VECTOR(Scalar, 2) Key;

  /**
   *  给定相机内外参数与世界座标系下的三维点座标，
   *  计算包括镜头畸变的投影后的影像座标，
   */
  static Key WorldPointProjectToImageKey(
    const IntrinsicParams& intrinsic_params,
    const ExtrinsicParams& extrinsic_params,
    const Point3D& point)
  {
    Point3D camera_point = point - extrinsic_params.position();
    camera_point = extrinsic_params.rotation() * camera_point;
    camera_point /= camera_point[2];

    Key radial_delta_key;
    RadialDistortor<Scalar>()(intrinsic_params.k1(),
                              intrinsic_params.k2(),
                              intrinsic_params.k3(),
                              camera_point[0],
                              camera_point[1],
                              radial_delta_key[0],
                              radial_delta_key[1]);

    Key decentering_delta_key;
    DecenteringDistortor<Scalar>()(intrinsic_params.d1(),
                                   intrinsic_params.d2(),
                                   camera_point[0],
                                   camera_point[1],
                                   decentering_delta_key[0],
                                   decentering_delta_key[1]);

    camera_point.segment(0, 2) += radial_delta_key;
    camera_point.segment(0, 2) += decentering_delta_key;
    camera_point = intrinsic_params.GetKMatrix() * camera_point;

    return camera_point.segment(0, 2);
  }

  /**
   *  给定相机内外参数与世界座标系下的三维点座标，
   *  计算不包括镜头畸变的投影后的影像座标。
   */
  static Key WorldPointProjectToImageKeyNoDistort(
    const IntrinsicParams& intrinsic_params,
    const ExtrinsicParams& extrinsic_params,
    const Point3D& point)
  {
    Point3D camera_point = point - extrinsic_params.position();
    camera_point = extrinsic_params.rotation() * camera_point;
    camera_point = intrinsic_params.GetKMatrix() * camera_point;
    camera_point /= camera_point[2];

    return camera_point.segment(0, 2);
  }

};

}
}

#endif
