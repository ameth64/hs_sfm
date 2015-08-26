set(HOMOGRAPHY_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_homography2d_analytical_jacobian_matrix_calculator.cpp"
  "test_homography2d_vector_function.cpp"
  "test_homography2d_dlt_calculator.cpp"
  )
hslib_add_utest(homography SOURCES ${HOMOGRAPHY_UTEST_SOURCES})
yong_add_dependence(hs_sfm homography_utest
  DEPENDENT_PROJECT boost
    HEADER_ONLY
  DEPENDENT_PROJECT cereal
    HEADER_ONLY
  DEPENDENT_PROJECT eigen
    HEADER_ONLY
    #DEPENDENT_PROJECT hs_fit
    #  HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    HEADER_ONLY
    #DEPENDENT_PROJECT hs_optimizor
    #  HEADER_ONLY
  )
hslib_add_mkl_dep(homography_utest)
hslib_add_utest_end(homography)
