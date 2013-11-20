#ifndef _HS_SFM_SFM_FILE_IO_EXTRINSIC_PARAMS_SET_LOADER_HPP_
#define _HS_SFM_SFM_FILE_IO_EXTRINSIC_PARAMS_SET_LOADER_HPP_

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
struct ExtrinsicParamsSetLoader
{
  typedef _Scalar Scalar;
  typedef int Err;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;

  Err operator() (const std::string& extrinsic_set_path,
                  ExtrinsicParamsContainer& extrinsic_params_set) const
  {
    std::ifstream extrinsic_set_file(extrinsic_set_path.c_str(),
                                     std::ios::in);
    if (!extrinsic_set_file.is_open())
    {
      return -1;
    }

    size_t number_of_cameras;
    extrinsic_set_file>>number_of_cameras;
    extrinsic_params_set.resize(number_of_cameras);
    for (size_t i = 0; i < number_of_cameras; i++)
    {
      extrinsic_set_file>>extrinsic_params_set[i].rotation()[0]
                        >>extrinsic_params_set[i].rotation()[1]
                        >>extrinsic_params_set[i].rotation()[2]
                        >>extrinsic_params_set[i].position()[0]
                        >>extrinsic_params_set[i].position()[1]
                        >>extrinsic_params_set[i].position()[2];
    }

    return 0;
  }
};

}
}
}

#endif
