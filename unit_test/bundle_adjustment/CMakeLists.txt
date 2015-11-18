set(BA_UTEST_SOURCES
  "${PROJECT_SOURCE_DIR}/unit_test/test_main.cpp"
  #"test_ba_naive_normal_mle_meta.hpp"
  #"test_ba_naive_jacobian_matrix_calculator.cpp"
  #"test_ba_naive_monte_carlo_simulate.cpp"
  #"test_ba_naive_normal_equation_builder.cpp"
  #"test_ba_naive_normal_equation_solver.cpp"
  #"test_ba_naive_vector_function.cpp"
  #"test_ba_naive_ground_truth.cpp"
  #"test_ba_gcp_constrained_gauge.cpp"
  #"test_ba_gcp_constrained_jacobian_matrix_calculator.cpp"
  #"test_ba_gcp_constrained_normal_equation_builder.cpp"
  #"test_ba_gcp_constrained_vector_function.cpp"
  "test_camera_shared_vector_function.cpp"
  "test_camera_shared_jacobian_matrix_calculator.cpp"
  "test_camera_shared_analytical_derivative.cpp"
  "test_camera_shared_normal_equation_builder.cpp"
  "test_camera_shared_normal_equation_solver.cpp"
  #"test_camera_shared_monte_carlo_simulate.cpp"
  "test_camera_shared_levenberg_marquardt_optimizor.cpp"
  "test_camera_shared_relative_pair.cpp"
  )
hslib_add_utest(bundle_adjustment SOURCES ${BA_UTEST_SOURCES})
if (MSVC)
  set_source_files_properties(
    "test_camera_shared_relative_pair.cpp"
    PROPERTIES COMPILE_FLAGS "/bigobj")
endif()
yong_add_dependence(hs_sfm bundle_adjustment_utest
  DEPENDENT_PROJECT boost
    HEADER_ONLY
  DEPENDENT_PROJECT cereal
    HEADER_ONLY
  DEPENDENT_PROJECT ceres
    DEPENDENT_LIBS ceres
  DEPENDENT_PROJECT eigen
    HEADER_ONLY
  DEPENDENT_PROJECT jpeg_turbo
    DEPENDENT_LIBS jpeg_turbo
  DEPENDENT_PROJECT png
    DEPENDENT_LIBS png
  DEPENDENT_PROJECT SuiteSparse
    DEPENDENT_LIBS amd camd cholmod colamd suitesparseconfig
  DEPENDENT_PROJECT tiff
    DEPENDENT_LIBS tiff
  DEPENDENT_PROJECT zlib
    DEPENDENT_LIBS zlib
  DEPENDENT_PROJECT hs_image_io
    DEPENDENT_LIBS whole_io
  DEPENDENT_PROJECT hs_math
    HEADER_ONLY
  DEPENDENT_PROJECT hs_optimizor
    HEADER_ONLY
  DEPENDENT_PROJECT hs_test_utility
    HEADER_ONLY
  )
hslib_add_mkl_dep(bundle_adjustment_utest)
hslib_add_utest_end(bundle_adjustment)
