cmake_minimum_required(VERSION 2.8)

include("cmake/hslib.cmake")
include("cmake/hsprj.cmake")

set(HSPRJ_NAME "hs_sfm")
hs_project(${HSPRJ_NAME})

#bundle adjustment
set(BA_NAME "bundle_adjustment")
set(BA_SOURCES
  )
hslib_add_library(${BA_NAME} ${HSPRJ_NAME} ${BA_SOURCES})

set(BA_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_base.hpp"
  "unit_test/bundle_adjustment/test_ba_naive_jac.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_nllso_normal_equation_builder.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_vec_func.cpp"
  )
hslib_add_utest(${BA_NAME} SOURCES ${BA_UTEST_SOURCES})
hslib_add_dep(${BA_NAME}_utest
  3RDPARTY_DEPS eigen HEADONLY
  3RDPARTY_DEPS boost HEADONLY
  HS_DEPS hs_math HEADONLY
  )
hslib_post_build(${BA_NAME}_utest)

#utility
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
