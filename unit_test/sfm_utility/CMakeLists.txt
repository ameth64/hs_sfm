set(UTILITY_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_camera_type.cpp"
  "test_key_type.cpp"
  "test_match_type.cpp"
  "test_ransac_fit_plane.cpp"
  "test_similar_transform_estimator.cpp"
  "test_matches_tracks_convertor.cpp"
  "test_undistortor.cpp"
  )
hslib_add_utest(sfm_utility SOURCES ${UTILITY_UTEST_SOURCES})
yong_add_dependence(hs_sfm sfm_utility_utest
  DEPENDENT_PROJECT boost
    COMMIT_HASH ${BOOST_HASH}
    HEADER_ONLY
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
  DEPENDENT_PROJECT hs_test_utility
    COMMIT_HASH ${HS_TEST_UTILITY_HASH}
    HEADER_ONLY
  )
hslib_add_mkl_dep(sfm_utility_utest)
hslib_add_utest_end(sfm_utility)
