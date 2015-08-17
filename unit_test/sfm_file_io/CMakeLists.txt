set(FILE_IO_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_loader_saver_consistency.cpp"
  )
hslib_add_utest(sfm_file_io SOURCES ${FILE_IO_UTEST_SOURCES})
yong_add_dependence(hs_sfm sfm_file_io_utest
  DEPENDENT_PROJECT boost
    COMMIT_HASH ${boost_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT cereal
    COMMIT_HASH ${cereal_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT eigen
    COMMIT_HASH ${eigen_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    COMMIT_HASH ${hs_math_HASH}
    HEADER_ONLY
    #DEPENDENT_PROJECT hs_test_utility
    #  COMMIT_HASH ${HS_TEST_UTILITY_HASH}
    #  HEADER_ONLY
  )
hslib_add_mkl_dep(sfm_file_io_utest)
hslib_add_utest_end(sfm_file_io)
