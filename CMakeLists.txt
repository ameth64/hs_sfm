cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0054 OLD)

include("cmake/yong_injector/injector.cmake")
include("cmake/yong_hsprj.cmake")
include("cmake/yong_hslib.cmake")

hs_project(hs_sfm)

yong_option(BUILD_UNIT_TEST OFF)

hs_parse_dependencies("dependencies.txt")

if (BUILD_UNIT_TEST)
  add_subdirectory("unit_test")
endif()

hs_project_end(hs_sfm)
