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
  "unit_test/bundle_adjustment/test_ba_naive_jacobian_matrix_calculator.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_monte_carlo_simulate.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_normal_equation_builder.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_normal_equation_solver.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_vector_function.cpp"
  "unit_test/bundle_adjustment/test_ba_naive_ground_truth.cpp"
  "unit_test/bundle_adjustment/test_ba_gcp_constrained_gauge.cpp"
  "unit_test/bundle_adjustment/test_ba_gcp_constrained_jacobian_matrix_calculator.cpp"
  "unit_test/bundle_adjustment/test_ba_gcp_constrained_normal_equation_builder.cpp"
  "unit_test/bundle_adjustment/test_ba_gcp_constrained_vector_function.cpp"
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
set(UTILITY_NAME "sfm_utility")
set(UTILITY_SOURCES
  )
hslib_add_library(${UTILITY_NAME} ${HSPRJ_NAME} ${UTILITY_SOURCES})

set(UTILITY_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/sfm_utility/test_camera_type.cpp"
  "unit_test/sfm_utility/test_ransac_fit_plane.cpp"
  "unit_test/sfm_utility/test_similar_transform_estimator.cpp"
  "unit_test/sfm_utility/test_matches_tracks_convertor.cpp"
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

#incremental
set(INCREMENTAL_NAME "incremental")
set(INCREMENTAL_SOURCES
  )
hslib_add_library(${INCREMENTAL_NAME} ${HSPRJ_NAME} ${INCREMENTAL_SOURCES})

set(INCREMENTAL_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/incremental/test_incremental.cpp"
  "unit_test/incremental/test_scene_expandor.cpp"
  "unit_test/incremental/test_bundle_adjustment_gcp_constrained_optimizor.cpp"
  #"unit_test/incremental/test_scene_expandor_gcp_constrained.cpp"
  )
hslib_add_utest(${INCREMENTAL_NAME} SOURCES ${INCREMENTAL_UTEST_SOURCES})
hslib_add_dep(${INCREMENTAL_NAME}_utest
  3RDPARTY_DEPS eigen HEADONLY
  HS_DEPS hs_math HEADONLY
  HS_DEPS hs_test_utility DEP_LIBS test_env)
hslib_post_build(${INCREMENTAL_NAME}_utest)

#essential
set(ESSENTIAL_NAME "essential")
set(ESSENTIAL_SOURCES
  )
hslib_add_library(${ESSENTIAL_NAME} ${HSPRJ_NAME} ${ESSENTIAL_SOURCES})

set(ESSENTIAL_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/essential/test_ematrix_5_points_calculator.cpp"
  "unit_test/essential/test_ematrix_5_points_ransac_refiner.cpp"
  "unit_test/essential/test_ematrix_extrinsic_params_points_calculator.cpp"
  )
hslib_add_utest(${ESSENTIAL_NAME} SOURCES ${ESSENTIAL_UTEST_SOURCES})
hslib_add_dep(${ESSENTIAL_NAME}_utest
  3RDPARTY_DEPS eigen HEADONLY
  HS_DEPS hs_math HEADONLY)
hslib_post_build(${ESSENTIAL_NAME}_utest)

#fundamental
set(FUNDAMENTAL_NAME "fundamental")
set(FUNDAMENTAL_SOURCES
  )
hslib_add_library(${FUNDAMENTAL_NAME} ${HSPRJ_NAME} ${FUNDAMENTAL_SOURCES})

set(FUNDAMENTAL_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/fundamental/test_linear_8_points_calculator.cpp"
  "unit_test/fundamental/test_linear_8_points_ransac_refiner.cpp"
  )
hslib_add_utest(${FUNDAMENTAL_NAME} SOURCES ${FUNDAMENTAL_UTEST_SOURCES})
hslib_add_dep(${FUNDAMENTAL_NAME}_utest
  3RDPARTY_DEPS eigen HEADONLY
  HS_DEPS hs_math HEADONLY)
hslib_post_build(${FUNDAMENTAL_NAME}_utest)

#homography
set(HOMOGRAPHY_NAME "homography")
set(HOMOGRAPHY_SOURCES
  )
hslib_add_library(${HOMOGRAPHY_NAME} ${HSPRJ_NAME} ${HOMOGRAPHY_SOURCES})

set(HOMOGRAPHY_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/homography/test_homography2d_analytical_jacobian_matrix_calculator.cpp"
  "unit_test/homography/test_homography2d_vector_function.cpp"
  )
hslib_add_utest(${HOMOGRAPHY_NAME} SOURCES ${HOMOGRAPHY_UTEST_SOURCES})
hslib_add_dep(${HOMOGRAPHY_NAME}_utest
  3RDPARTY_DEPS eigen HEADONLY
  HS_DEPS hs_math HEADONLY)
hslib_post_build(${HOMOGRAPHY_NAME}_utest)

#projective
set(PROJECTIVE_NAME "projective")
set(PROJECTIVE_SOURCES
  )
hslib_add_library(${PROJECTIVE_NAME} ${HSPRJ_NAME} ${PROJECTIVE_SOURCES})

set(PROJECTIVE_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/projective/test_pmatrix_dlt_calculator.cpp"
  "unit_test/projective/test_pmatrix_dlt_ransac_refiner.cpp"
  "unit_test/projective/test_single_camera_params_maximum_likelihood_estimator.cpp"
  "unit_test/projective/mle/test_intrinsic_constrained_analytical_jacobian_matrix_calculator.cpp"
  "unit_test/projective/mle/test_intrinsic_constrained_vector_function.cpp"
  "unit_test/projective/mle/test_intrinsic_constrained_normal_equation_builder.cpp"
  )
hslib_add_utest(${PROJECTIVE_NAME} SOURCES ${PROJECTIVE_UTEST_SOURCES})
hslib_add_dep(${PROJECTIVE_NAME}_utest
  3RDPARTY_DEPS eigen HEADONLY
  HS_DEPS hs_math HEADONLY
  HS_DEPS hs_optimizor HEADONLY)
hslib_post_build(${PROJECTIVE_NAME}_utest)

#synthetic
set(SYNTHETIC_NAME "synthetic")
set(SYNTHETIC_SOURCES
  )
hslib_add_library(${SYNTHETIC_NAME} ${HSPRJ_NAME} ${SYNTHETIC_SOURCES})

set(SYNTHETIC_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/synthetic/test_keyset_generator.cpp"
  "unit_test/synthetic/test_relative_generator.cpp"
  "unit_test/synthetic/test_scene_generator.cpp"
  )
hslib_add_utest(${SYNTHETIC_NAME} SOURCES ${SYNTHETIC_UTEST_SOURCES})
hslib_add_dep(${SYNTHETIC_NAME}_utest
  3RDPARTY_DEPS eigen HEADONLY
  HS_DEPS hs_math HEADONLY)
hslib_post_build(${SYNTHETIC_NAME}_utest)

#sfm file io
set(FILEIO_NAME "sfm_file_io")
set(FILEIO_SOUCES 
  )
hslib_add_library(${FILEIO_NAME} ${HSPRJ_NAME} ${FILEIO_SOURCES})

set(FILEIO_UTEST_SOURCES
  "unit_test/test_main.cpp"
  "unit_test/sfm_file_io/test_loader_saver_consistency.cpp"
  )
hslib_add_utest(${FILEIO_NAME} SOURCES ${FILEIO_UTEST_SOURCES})
hslib_add_dep(${FILEIO_NAME}_utest
  3RDPARTY_DEPS eigen HEADONLY
  3RDPARTY_DEPS boost HEADONLY
  HS_DEPS hs_math HEADONLY
  HS_DEPS hs_test_utility DEP_LIBS test_env
  )
hslib_post_build(${FILEIO_NAME}_utest)
