set(ESSENTIAL_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_ematrix_5_points_calculator.cpp"
  "test_ematrix_5_points_ransac_refiner.cpp"
  "test_ematrix_extrinsic_params_points_calculator.cpp"
  )
hslib_add_utest(essential SOURCES ${ESSENTIAL_UTEST_SOURCES})
yong_add_dependence(hs_sfm essential_utest
  DEPENDENT_PROJECT boost
    HEADER_ONLY
  DEPENDENT_PROJECT cereal
    HEADER_ONLY
  DEPENDENT_PROJECT eigen
    HEADER_ONLY
  DEPENDENT_PROJECT hs_fit
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    HEADER_ONLY
  DEPENDENT_PROJECT hs_optimizor
    HEADER_ONLY
  )
hslib_add_mkl_dep(essential_utest)
hslib_add_utest_end(essential)
