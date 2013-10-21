#include <gtest/gtest.h>

#include "hs_sfm/sfm_utility/synthetic_scene_generator.hpp"

#include "hs_sfm/fundamental/linear_8_points_calculator.hpp"

namespace
{

template <typename _Scalar, typename _ImageDimension>
class TestLinear8PointsCalculator
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
};

}