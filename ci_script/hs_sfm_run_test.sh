#!/bin/sh

. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/initialize.sh "$@"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/run_test.sh \
  "Release" "unit_test/bundle_adjustment" "bundle_adjustment_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/run_test.sh \
  "Release" "unit_test/calibrate" "calibrate_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/run_test.sh \
  "Release" "unit_test/essential" "essential_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/run_test.sh \
  "Release" "unit_test/fundamental" "fundamental_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/run_test.sh \
  "Release" "unit_test/homography" "homography_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/run_test.sh \
  "Release" "unit_test/projective" "projective_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/run_test.sh \
  "Release" "unit_test/sfm_file_io" "sfm_file_io_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/run_test.sh \
  "Release" "unit_test/sfm_pipeline" "sfm_pipeline_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/run_test.sh \
  "Release" "unit_test/sfm_utility" "sfm_utility_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/run_test.sh \
  "Release" "unit_test/synthetic" "synthetic_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/run_test.sh \
  "Release" "unit_test/triangulate" "_triangulateutest"
