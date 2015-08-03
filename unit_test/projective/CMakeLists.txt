set(PROJECTIVE_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_pmatrix_dlt_calculator.cpp"
  "test_pmatrix_dlt_ransac_refiner.cpp"
  "test_single_camera_params_maximum_likelihood_estimator.cpp"
  "mle/test_intrinsic_constrained_analytical_jacobian_matrix_calculator.cpp"
  "mle/test_intrinsic_constrained_vector_function.cpp"
  "mle/test_intrinsic_constrained_normal_equation_builder.cpp"
  )
hslib_add_utest(projective SOURCES ${PROJECTIVE_UTEST_SOURCES})
yong_add_dependence(hs_sfm projective_utest
  DEPENDENT_PROJECT cereal
    COMMIT_HASH ${CEREAL_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT ceres
    COMMIT_HASH ${CERES_HASH}
    DEPENDENT_LIBS ceres
  DEPENDENT_PROJECT eigen
    COMMIT_HASH ${EIGEN_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_fit
    COMMIT_HASH ${HS_FIT_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    COMMIT_HASH ${HS_MATH_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_optimizor
    COMMIT_HASH ${HS_OPTIMIZOR_HASH}
    HEADER_ONLY
  )
hslib_add_mkl_dep(projective_utest)
hslib_add_utest_end(projective)
