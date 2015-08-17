#!/bin/sh

. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/initialize.sh "$@"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "Release" "bundle_adjustment_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "Release" "calibrate_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "Release" "essential_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "Release" "fundamental_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "Release" "homography_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "Release" "projective_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "Release" "sfm_file_io_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "Release" "sfm_pipeline_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "Release" "sfm_utility_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "Release" "synthetic_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "Release" "triangulate_utest"
