#ifndef _HS_SFM_SFM_FILE_IO_INTRINSIC_PARAMS_SET_SAVER_HPP_
#define _HS_SFM_SFM_FILE_IO_INTRINSIC_PARAMS_SET_SAVER_HPP_

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
struct IntrinsicParamsSetSaver
{
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;

  Err operator() (const std::string& intrinsic_set_path,
                  const IntrinsicParamsContainer& intrinsic_params_set) const
  {
    std::ofstream intrinsic_set_file(intrinsic_set_path.c_str(),
                                     std::ios::out);
    if (!intrinsic_set_file.is_open())
    {
      return -1;
    }

    size_t number_of_cameras = intrinsic_params_set.size();
    intrinsic_set_file<<number_of_cameras<<"\n";
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      const IntrinsicParams& intrinsic_params = intrinsic_params_set[i];
      intrinsic_set_file<<intrinsic_params.focal_length()<<"\n";
      intrinsic_set_file<<intrinsic_params.skew()<<"\n";
      intrinsic_set_file<<intrinsic_params.principal_point_x()<<"\n";
      intrinsic_set_file<<intrinsic_params.principal_point_y()<<"\n";
      intrinsic_set_file<<intrinsic_params.pixel_ratio()<<"\n";
    }

    return 0;
  }
};

}
}
}

#endif
