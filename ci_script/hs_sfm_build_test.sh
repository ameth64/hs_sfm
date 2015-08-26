#!/bin/bash

. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/initialize.sh "$@"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "bundle_adjustment_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "calibrate_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "essential_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "fundamental_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "homography_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "projective_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "sfm_file_io_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "sfm_pipeline_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "sfm_utility_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "synthetic_utest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/build_test.sh \
  "triangulate_utest"
