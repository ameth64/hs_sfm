set(CALIBRATE_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_planar_calibrator.cpp"
  )
hslib_add_utest(calibrate SOURCES ${CALIBRATE_UTEST_SOURCES})
yong_add_dependence(hs_sfm calibrate_utest
  DEPENDENT_PROJECT boost
    HEADER_ONLY
  DEPENDENT_PROJECT cereal
    HEADER_ONLY
  DEPENDENT_PROJECT eigen
    HEADER_ONLY
  DEPENDENT_PROJECT OpenCV
    DEPENDENT_LIBS opencv_core
                   opencv_calib3d
                   opencv_imgproc
                   opencv_features2d
                   opencv_highgui
                   opencv_flann
  DEPENDENT_PROJECT zlib
    DEPENDENT_LIBS zlib
  DEPENDENT_PROJECT hs_math
    HEADER_ONLY
  )
hslib_add_mkl_dep(calibrate_utest)
hslib_add_utest_end(calibrate)
