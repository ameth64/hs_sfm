#ifndef _HS_SFM_SFM_FILE_IO_POINTS_SAVER_HPP_
#define _HS_SFM_SFM_FILE_IO_POINTS_SAVER_HPP_

#include <string>
#include <fstream>

#include "hs_math/linear_algebra/eigen_macro.hpp"

namespace hs
{
namespace sfm
{
namespace fileio
{

template <typename _Scalar>
struct PointsSaver
{
  typedef _Scalar Scalar;
  typedef int Err;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

  Err operator() (const std::string& points_path,
                  const PointContainer& points) const
  {
    std::ofstream points_file(points_path.c_str(), std::ios::out);
    if (!points_file.is_open())
    {
      return -1;
    }

    size_t number_of_points = points.size();
    points_file<<number_of_points<<"\n";
    for (size_t i = 0; i < number_of_points; i++)
    {
      points_file<<points[i][0]<<" "
                 <<points[i][1]<<" "
                 <<points[i][2]<<"\n";
    }

    return 0;
  }
};

}
}
}

#endif
