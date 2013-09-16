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
  "unit_test/bundle_adjustment/test_ba_naive_normal_mle_meta.hpp"
  "unit_test/bundle_adjustment/test_ba_naive_jac.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_monte_carlo_simulate.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_nllso_levenberg_marquardt.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_nllso_normal_equation_builder.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_nllso_normal_equation_solver.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_vec_func.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_xcov_calculator.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_ground_truth.cpp"
  )
hslib_add_utest(${BA_NAME} SOURCES ${BA_UTEST_SOURCES})
hslib_add_dep(${BA_NAME}_utest
  3RDPARTY_DEPS eigen HEADONLY
  3RDPARTY_DEPS boost HEADONLY
  HS_DEPS hs_math HEADONLY
  HS_DEPS hs_optimizor HEADONLY
  )
hslib_post_build(${BA_NAME}_utest)

#utility
set(UTILITY_NAME "utility")
set(UTILITY_SOURCES
  )
hslib_add_library(${UTILITY_NAME} ${HSPRJ_NAME} ${UTILITY_SOURCES})

set(UTILITY_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/utility/test_camera_type.cpp"
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

#triangulate
set(TRIANGULATE_NAME "triangulate")
set(TRIANGULATE_SOURCES
  )
hslib_add_library(${TRIANGULATE_NAME} ${HSPRJ_NAME} ${TRIANGULATE_SOURCES})

set(TRIANGULATE_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/triangulate/test_multiple_view_base.hpp"
  "unit_test/triangulate/test_multiple_view_normal_mle_meta.hpp"
  "unit_test/triangulate/test_multiple_view_vector_function.cpp"
  "unit_test/triangulate/test_multiple_view_monte_carlo_simulate.cpp"
  )
hslib_add_utest(${TRIANGULATE_NAME} SOURCES ${TRIANGULATE_UTEST_SOURCES})
hslib_add_dep(${TRIANGULATE_NAME}_utest
  3RDPARTY_DEPS eigen HEADONLY
  HS_DEPS hs_math HEADONLY)
hslib_post_build(${TRIANGULATE_NAME}_utest)
