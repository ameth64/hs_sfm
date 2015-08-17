set(FUNDAMENTAL_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_linear_8_points_calculator.cpp"
  "test_linear_8_points_ransac_refiner.cpp"
  )
hslib_add_utest(fundamental SOURCES ${FUNDAMENTAL_UTEST_SOURCES})
yong_add_dependence(hs_sfm fundamental_utest
  DEPENDENT_PROJECT cereal
    COMMIT_HASH ${cereal_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT eigen
    COMMIT_HASH ${eigen_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_fit
    COMMIT_HASH ${hs_fit_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    COMMIT_HASH ${hs_math_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_optimizor
    COMMIT_HASH ${hs_optimizor_HASH}
    HEADER_ONLY
  )
hslib_add_mkl_dep(fundamental_utest)
hslib_add_utest_end(fundamental)
