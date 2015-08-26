#!/bin/bash

. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/initialize.sh "$@"

. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "yong" "boost"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "yong" "cereal"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "yong" "eigen"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "eigen" "yong" "ceres"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "yong" "flann"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "yong" "gtest"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" \
  --extra_cmake_defs "-DNASM=${CI_NASM_PATH}" "yong" "jpeg_turbo"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "yong" "lemon"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "yong" "osi_clp"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "yong" "zlib"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "zlib" "yong" "png"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "jpeg_turbo zlib" "yong" "tiff"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "jpeg_turbo png tiff zlib" "yong" "OpenCV"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "yong_injector" --deps "cereal ceres eigen flann jpeg_turbo lemon osi_clp png tiff zlib" \
  "yong" "openmvg"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "cmake ci_script/hslib_ci_script" \
  --deps "jpeg_turbo png tiff zlib gtest" "yong" "hs_image_io"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "boost cereal eigen gtest" "yong" "hs_math"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "boost cereal eigen gtest hs_math" "yong" "hs_fit"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "gtest" "yong" "hs_progress"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "boost cereal eigen gtest hs_math" \
  "yong" "hs_test_utility"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "boost cereal eigen gtest hs_math hs_test_utility" \
  "yong" "hs_optimizor"

. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/cmake_generate.sh \
  --extra_cmake_defs "-DHSLIB_BUILD_DOC=0 -DHSLIB_COPY_TEST_DATA=1 -DBUILD_UNIT_TEST=1" \
  --deps "boost eigen cereal ceres flann gtest jpeg_turbo lemon osi_clp zlib png tiff OpenCV openmvg hs_image_io hs_math hs_fit hs_progress hs_test_utility hs_optimizor"
