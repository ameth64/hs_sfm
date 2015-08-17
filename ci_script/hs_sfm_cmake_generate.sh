#!/bin/sh

. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/initialize.sh "$@"

. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "Release" "boost"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "Release" "cereal"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "Release" "eigen"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "eigen" "Release" "ceres"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "Release" "flann"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "Release" "gtest"
jpeg_turbo_extra_cmake_defs="-DNASM=nasm"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" \
  --extra_cmake_defs "${jpeg_turbo_extra_cmake_defs}" "Release" "jpeg_turbo"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "Release" "lemon"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "Release" "osi_clp"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" "Release" "zlib"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "zlib" "Release" "png"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "jpeg_turbo zlib" "Release" "tiff"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "jpeg_turbo png tiff zlib" "Release" "OpenCV"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "cereal ceres eigen flann jpeg_turbo lemon osi_clp png tiff zlib" \
  "Release" "openmvg"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "jpeg_turbo png tiff zlib gtest" "Release" "hs_image_io"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "boost cereal eigen gtest" "Release" "hs_math"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "boost cereal eigen gtest hs_math" "Release" "hs_fit"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "gtest" "Release" "hs_progress"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "boost cereal eigen gtest hs_math" \
  "Release" "hs_test_utility"
. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/resolve_dependency.sh \
  --submodules "ALL" --deps "boost cereal eigen gtest hs_math hs_test_utility" \
  "Release" "hs_optimizor"

. ${CI_PROJECT_DIR}/ci_script/hslib_ci_script/cmake_generate.sh
