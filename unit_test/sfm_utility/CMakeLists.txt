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
    COMMIT_HASH ${boost_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT cereal
    COMMIT_HASH ${cereal_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT ceres
    COMMIT_HASH ${ceres_HASH}
    DEPENDENT_LIBS ceres
  DEPENDENT_PROJECT eigen
    COMMIT_HASH ${eigen_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_fit
    COMMIT_HASH ${hs_fit_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    COMMIT_HASH ${hs_math_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_test_utility
    COMMIT_HASH ${hs_test_utility_HASH}
    HEADER_ONLY
  )
hslib_add_mkl_dep(sfm_utility_utest)
hslib_add_utest_end(sfm_utility)
