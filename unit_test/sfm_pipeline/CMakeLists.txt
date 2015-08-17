set(PIPELINE_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_incremental_sfm.cpp"
  "test_global_sfm.cpp"
  "test_point_cloud_norm_calculator.cpp"
  )
hslib_add_utest(sfm_pipeline SOURCES ${PIPELINE_UTEST_SOURCES})
yong_add_dependence(hs_sfm sfm_pipeline_utest
  DEPENDENT_PROJECT boost
    COMMIT_HASH ${boost_HASH}
    DEPENDENT_LIBS boost_filesystem boost_system
  DEPENDENT_PROJECT cereal
    COMMIT_HASH ${cereal_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT ceres
    COMMIT_HASH ${ceres_HASH}
    DEPENDENT_LIBS ceres
  DEPENDENT_PROJECT eigen
    COMMIT_HASH ${eigen_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT flann
    COMMIT_HASH ${flann_HASH}
    DEPENDENT_LIBS flann
  DEPENDENT_PROJECT lemon
    COMMIT_HASH ${lemon_HASH}
    DEPENDENT_LIBS lemon
  DEPENDENT_PROJECT openmvg
    COMMIT_HASH ${openmvg_HASH}
    DEPENDENT_LIBS openMVG_sfm
                   openMVG_multiview
                   openMVG_lInftyComputerVision
                   openMVG_numeric
                   openMVG_kvld
                   openMVG_features
                   openMVG_system
                   openMVG_stlplus
  DEPENDENT_PROJECT osi_clp
    COMMIT_HASH ${osi_clp_HASH}
    DEPENDENT_LIBS clp CoinUtils Osi OsiClpSolver
  DEPENDENT_PROJECT zlib
    COMMIT_HASH ${zlib_HASH}
    DEPENDENT_LIBS zlib
  DEPENDENT_PROJECT hs_fit
    COMMIT_HASH ${hs_fit_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    COMMIT_HASH ${hs_math_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_optimizor
    COMMIT_HASH ${hs_optimizor_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_progress
    COMMIT_HASH ${hs_progress_HASH}
    DEPENDENT_LIBS progress_utility
  DEPENDENT_PROJECT hs_test_utility
    COMMIT_HASH ${hs_test_utility_HASH}
    HEADER_ONLY
  )
hslib_add_mkl_dep(sfm_pipeline_utest)
hslib_add_utest_end(sfm_pipeline)
