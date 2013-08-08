cmake_minimum_required(VERSION 2.8)

include("cmake/hslib.cmake")
include("cmake/hsprj.cmake")

set(HSPRJ_NAME "hs_sfm")
hs_project(${HSPRJ_NAME})

set(UTILITY_NAME "utility")
set(UTILITY_SOURCES
  )
hslib_add_library(${UTILITY_NAME} ${HSPRJ_NAME} ${UTILITY_SOURCES})

set(UTILITY_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/utility/test_cam_type.cpp"
  "unit_test/utility/test_ransac_fit_plane.cpp"
  "unit_test/utility/test_sfm_file_io.cpp"
  "unit_test/utility/test_synthetic_scene_generator.cpp"
  )
hslib_add_utest(${UTILITY_NAME} SOURCES ${UTILITY_UTEST_SOURCES})
hslib_add_dep(${UTILITY_NAME}_utest
  3RDPARTY_DEPS eigen HEADONLY
  3RDPARTY_DEPS boost HEADONLY
  HS_DEPS hs_math HEADONLY
  HS_DEPS hs_test_utility DEP_LIBS test_env
  HS_DEPS hs_fit HEADONLY
  )
hslib_post_build(${UTILITY_NAME}_utest)
