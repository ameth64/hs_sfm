set(UTILITY_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_camera_type.cpp"
  "test_key_type.cpp"
  "test_match_type.cpp"
  "test_ransac_fit_plane.cpp"
  "test_similar_transform_estimator.cpp"
  "test_matches_tracks_convertor.cpp"
  "test_undistortor.cpp"
  "test_camera_rotated_bounding_box_intersector.cpp"
  )
hslib_add_utest(sfm_utility SOURCES ${UTILITY_UTEST_SOURCES})
yong_add_dependence(hs_sfm sfm_utility_utest
  DEPENDENT_PROJECT boost
    HEADER_ONLY
  DEPENDENT_PROJECT cereal
    HEADER_ONLY
  DEPENDENT_PROJECT ceres
    DEPENDENT_LIBS ceres
  DEPENDENT_PROJECT eigen
    HEADER_ONLY
  DEPENDENT_PROJECT hs_fit
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    HEADER_ONLY
  DEPENDENT_PROJECT hs_test_utility
    HEADER_ONLY
  )
hslib_add_mkl_dep(sfm_utility_utest)
hslib_add_utest_end(sfm_utility)
