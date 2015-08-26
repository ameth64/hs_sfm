set(FUNDAMENTAL_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_linear_8_points_calculator.cpp"
  "test_linear_8_points_ransac_refiner.cpp"
  )
hslib_add_utest(fundamental SOURCES ${FUNDAMENTAL_UTEST_SOURCES})
yong_add_dependence(hs_sfm fundamental_utest
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
hslib_add_mkl_dep(fundamental_utest)
hslib_add_utest_end(fundamental)
