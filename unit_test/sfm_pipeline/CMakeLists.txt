set(PIPELINE_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_incremental_sfm.cpp"
  #"test_global_sfm.cpp"
  "test_point_cloud_norm_calculator.cpp"
  )
if (MSVC)
  set_source_files_properties(
    "test_incremental_sfm.cpp"
    PROPERTIES COMPILE_FLAGS "/bigobj")
endif()
hslib_add_utest(sfm_pipeline SOURCES ${PIPELINE_UTEST_SOURCES})
yong_add_dependence(hs_sfm sfm_pipeline_utest
  DEPENDENT_PROJECT boost
    DEPENDENT_LIBS boost_filesystem boost_system
  DEPENDENT_PROJECT cereal
    HEADER_ONLY
  DEPENDENT_PROJECT ceres
    DEPENDENT_LIBS ceres
  DEPENDENT_PROJECT eigen
    HEADER_ONLY
  DEPENDENT_PROJECT flann
    DEPENDENT_LIBS flann
  DEPENDENT_PROJECT lemon
    DEPENDENT_LIBS lemon
  DEPENDENT_PROJECT openmvg
    DEPENDENT_LIBS openMVG_sfm
                   openMVG_multiview
                   openMVG_lInftyComputerVision
                   openMVG_numeric
                   openMVG_kvld
                   openMVG_features
                   openMVG_system
                   openMVG_stlplus
  DEPENDENT_PROJECT osi_clp
    DEPENDENT_LIBS clp CoinUtils Osi OsiClpSolver
  DEPENDENT_PROJECT SuiteSparse
    DEPENDENT_LIBS amd camd cholmod colamd suitesparseconfig
  DEPENDENT_PROJECT zlib
    DEPENDENT_LIBS zlib
  DEPENDENT_PROJECT hs_fit
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    HEADER_ONLY
  DEPENDENT_PROJECT hs_optimizor
    HEADER_ONLY
  DEPENDENT_PROJECT hs_progress
    DEPENDENT_LIBS progress_utility
  DEPENDENT_PROJECT hs_test_utility
    HEADER_ONLY
  )
hslib_add_mkl_dep(sfm_pipeline_utest)
hslib_add_utest_end(sfm_pipeline)
