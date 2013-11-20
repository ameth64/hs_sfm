#ifndef _HS_SFM_SFM_FILE_IO_EXTRINSIC_PARAMS_SET_SAVER_HPP_
#define _HS_SFM_SFM_FILE_IO_EXTRINSIC_PARAMS_SET_SAVER_HPP_

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
struct ExtrinsicParamsSetSaver
{
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;

  Err operator() (const std::string& extrinsic_set_path,
                  const ExtrinsicParamsContainer& extrinsic_params_set) const
  {
    std::ofstream extrinsic_set_file(extrinsic_set_path.c_str(),
                                     std::ios::out);
    if (!extrinsic_set_file.is_open())
    {
      return -1;
    }

    size_t number_of_cameras = extrinsic_params_set.size();
    extrinsic_set_file<<number_of_cameras<<"\n";
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      extrinsic_set_file<<extrinsic_params_set[i].rotation()[0]<<" "
                        <<extrinsic_params_set[i].rotation()[1]<<" "
                        <<extrinsic_params_set[i].rotation()[2]<<" "
                        <<extrinsic_params_set[i].position()[0]<<" "
                        <<extrinsic_params_set[i].position()[1]<<" "
                        <<extrinsic_params_set[i].position()[2]<<"\n";
    }

    return 0;
  }
};

}
}
}

#endif
