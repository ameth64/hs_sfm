cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0054 OLD)

include("cmake/yong_injector/injector.cmake")
include("cmake/yong_hsprj.cmake")
include("cmake/yong_hslib.cmake")

hs_project(hs_sfm)

hs_parse_dependencies("dependencies.txt")

add_subdirectory("unit_test")

hs_project_end(hs_sfm)
