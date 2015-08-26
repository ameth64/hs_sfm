set(SYNTHETIC_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_flight_generator.cpp"
  "test_multiple_camera_generator.cpp"
  "test_keyset_generator.cpp"
  "test_relative_generator.cpp"
  "test_scene_generator.cpp"
  )
hslib_add_utest(synthetic SOURCES ${SYNTHETIC_UTEST_SOURCES})
yong_add_dependence(hs_sfm synthetic_utest
  DEPENDENT_PROJECT boost
    HEADER_ONLY
  DEPENDENT_PROJECT cereal
    HEADER_ONLY
  DEPENDENT_PROJECT eigen
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    HEADER_ONLY
  )
hslib_add_mkl_dep(synthetic_utest)
hslib_add_utest_end(synthetic)
