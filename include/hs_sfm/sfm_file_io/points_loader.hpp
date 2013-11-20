#ifndef _HS_SFM_SFM_FILE_IO_POINTS_LOADER_HPP_
#define _HS_SFM_SFM_FILE_IO_POINTS_LOADER_HPP_

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
struct PointsLoader
{
  typedef _Scalar Scalar;
  typedef int Err;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;

  Err operator() (const std::string& points_path,
                  PointContainer& points) const
  {
    std::ifstream points_file(points_path.c_str(), std::ios::in);
    if (!points_file.is_open())
    {
      return -1;
    }

    size_t number_of_points;
    points_file>>number_of_points;
    points.resize(number_of_points);
    for (size_t i = 0; i < number_of_points; i++)
    {
      points_file>>points[i][0]>>points[i][1]>>points[i][2];
    }

    return 0;
  }
};

}
}
}

#endif
