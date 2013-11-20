#ifndef _HS_SFM_SFM_FILE_IO_INTRINSIC_PARAMS_SET_LOADER_HPP_
#define _HS_SFM_SFM_FILE_IO_INTRINSIC_PARAMS_SET_LOADER_HPP_

#include <string>
#include <fstream>

#include "hs_sfm/sfm_utility/camera_type.hpp"

namespace hs
{
namespace sfm
{
namespace fileio
{

template <typename _Scalar>
struct IntrinsicParamsSetLoader
{
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;

  Err operator() (const std::string& intrinsic_set_path,
                  IntrinsicParamsContainer& intrinsic_params_set) const
  {
    std::ifstream intrinsic_set_file(intrinsic_set_path.c_str(),
                                     std::ios::in);
    if (!intrinsic_set_file.is_open())
    {
      return -1;
    }

    size_t number_of_cameras;
    intrinsic_set_file>>number_of_cameras;
    intrinsic_params_set.resize(number_of_cameras);
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      Scalar focal_length, skew,
             principal_point_x, principal_point_y,
             pixel_ratio;

      intrinsic_set_file>>focal_length
                        >>skew
                        >>principal_point_x
                        >>principal_point_y
                        >>pixel_ratio;

      intrinsic_params_set[i].set_focal_length(focal_length);
      intrinsic_params_set[i].set_skew(skew);
      intrinsic_params_set[i].set_principal_point_x(principal_point_x);
      intrinsic_params_set[i].set_principal_point_y(principal_point_y);
      intrinsic_params_set[i].set_pixel_ratio(pixel_ratio);
    }

    return 0;
  }

};

}
}
}

#endif
