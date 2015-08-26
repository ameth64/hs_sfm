set(TRIANGULATE_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_multiple_view_base.hpp"
  "test_multiple_view_normal_mle_meta.hpp"
  "test_multiple_view_vector_function.cpp"
  "test_multiple_view_monte_carlo_simulate.cpp"
  #"test_multiple_view_maximum_likelihood_estimator.cpp"
  )
hslib_add_utest(triangulate SOURCES ${TRIANGULATE_UTEST_SOURCES})
yong_add_dependence(hs_sfm triangulate_utest
  DEPENDENT_PROJECT boost
    HEADER_ONLY
  DEPENDENT_PROJECT cereal
    HEADER_ONLY
  DEPENDENT_PROJECT eigen
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    HEADER_ONLY
  DEPENDENT_PROJECT hs_optimizor
    HEADER_ONLY
  DEPENDENT_PROJECT hs_test_utility
    HEADER_ONLY
  )
hslib_add_mkl_dep(triangulate_utest)
hslib_add_utest_end(triangulate)
