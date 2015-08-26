set(FILE_IO_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_loader_saver_consistency.cpp"
  )
hslib_add_utest(sfm_file_io SOURCES ${FILE_IO_UTEST_SOURCES})
yong_add_dependence(hs_sfm sfm_file_io_utest
  DEPENDENT_PROJECT boost
    HEADER_ONLY
  DEPENDENT_PROJECT cereal
    HEADER_ONLY
  DEPENDENT_PROJECT eigen
    HEADER_ONLY
  DEPENDENT_PROJECT hs_math
    HEADER_ONLY
    #DEPENDENT_PROJECT hs_test_utility
    #  HEADER_ONLY
  )
hslib_add_mkl_dep(sfm_file_io_utest)
hslib_add_utest_end(sfm_file_io)
