set(PIPELINE_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_incremental_sfm.cpp"
  "test_global_sfm.cpp"
  "test_point_cloud_norm_calculator.cpp"
  )
hslib_add_utest(sfm_pipeline SOURCES ${PIPELINE_UTEST_SOURCES})
yong_add_dependence(hs_sfm sfm_pipeline_utest
  DEPENDENT_PROJECT boost
    COMMIT_HASH ${BOOST_HASH}
    DEPENDENT_LIBS boost_filesystem boost_system
  DEPENDENT_PROJECT cereal
    COMMIT_HASH ${CEREAL_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT ceres
    COMMIT_HASH ${CERES_HASH}
    DEPENDENT_LIBS ceres
  DEPENDENT_PROJECT eigen
    COMMIT_HASH ${EIGEN_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT flann
    COMMIT_HASH ${FLANN_HASH}
    DEPENDENT_LIBS flann
  DEPENDENT_PROJECT lemon
    COMMIT_HASH ${LEMON_HASH}
    DEPENDENT_LIBS lemon
  DEPENDENT_PROJECT openmvg
    COMMIT_HASH ${OPENMVG_HASH}
    DEPENDENT_LIBS openMVG_sfm
                   openMVG_multiview
                   openMVG_lInftyComputerVision
                   openMVG_numeric
                   openMVG_kvld
                   openMVG_features
                   openMVG_system
                   openMVG_stlplus
  DEPENDENT_PROJECT osi_clp
    COMMIT_HASH ${OSI_CLP_HASH}
    DEPENDENT_LIBS clp CoinUtils Osi OsiClpSolver
  DEPENDENT_PROJECT zlib
    COMMIT_HASH ${ZLIB_HASH}
    DEPENDENT_LIBS zlib
  DEPENDENT_PROJECT hs_fit
    COMMIT_HASH ${HS_FIT_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    COMMIT_HASH ${HS_MATH_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_optimizor
    COMMIT_HASH ${HS_OPTIMIZOR_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_progress
    COMMIT_HASH ${HS_PROGRESS_HASH}
    DEPENDENT_LIBS progress_utility
  DEPENDENT_PROJECT hs_test_utility
    COMMIT_HASH ${HS_TEST_UTILITY_HASH}
    HEADER_ONLY
  )
hslib_add_mkl_dep(sfm_pipeline_utest)
hslib_add_utest_end(sfm_pipeline)
