set(HOMOGRAPHY_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_homography2d_analytical_jacobian_matrix_calculator.cpp"
  "test_homography2d_vector_function.cpp"
  "test_homography2d_dlt_calculator.cpp"
  )
hslib_add_utest(homography SOURCES ${HOMOGRAPHY_UTEST_SOURCES})
yong_add_dependence(hs_sfm homography_utest
  DEPENDENT_PROJECT cereal
    COMMIT_HASH ${cereal_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT eigen
    COMMIT_HASH ${eigen_HASH}
    HEADER_ONLY
    #DEPENDENT_PROJECT hs_fit
    #  COMMIT_HASH ${HS_FIT_HASH}
    #  HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    COMMIT_HASH ${hs_math_HASH}
    HEADER_ONLY
    #DEPENDENT_PROJECT hs_optimizor
    #  COMMIT_HASH ${HS_OPTIMIZOR_HASH}
    #  HEADER_ONLY
  )
hslib_add_mkl_dep(homography_utest)
hslib_add_utest_end(homography)
