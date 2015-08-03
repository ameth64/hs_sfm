set(CALIBRATE_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  "test_planar_calibrator.cpp"
  )
hslib_add_utest(calibrate SOURCES ${CALIBRATE_UTEST_SOURCES})
yong_add_dependence(hs_sfm calibrate_utest
  DEPENDENT_PROJECT boost
    COMMIT_HASH ${BOOST_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT cereal
    COMMIT_HASH ${CEREAL_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT eigen
    COMMIT_HASH ${EIGEN_HASH}
    HEADER_ONLY
  DEPENDENT_PROJECT OpenCV
    COMMIT_HASH ${OPENCV_HASH}
    DEPENDENT_LIBS opencv_core
                   opencv_calib3d
                   opencv_imgproc
                   opencv_features2d
                   opencv_highgui
                   opencv_flann
  DEPENDENT_PROJECT zlib
    COMMIT_HASH ${ZLIB_HASH}
    DEPENDENT_LIBS zlib
  DEPENDENT_PROJECT hs_math
    COMMIT_HASH ${HS_MATH_HASH}
    HEADER_ONLY
  )
hslib_add_mkl_dep(calibrate_utest)
hslib_add_utest_end(calibrate)
