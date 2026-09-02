add_executable(lens_tests test_smoke.cpp test_image.cpp test_transfer.cpp test_cie.cpp test_upsample.cpp test_spectable.cpp test_dispersion.cpp test_warp.cpp test_vignette.cpp test_zernike.cpp test_pupil.cpp test_fft.cpp test_psf.cpp test_eff.cpp test_psfrings.cpp test_targets.cpp test_metrics.cpp test_pipeline.cpp test_acceptance.cpp test_lensfile.cpp test_colorfidelity.cpp test_flatfield.cpp test_grain.cpp)
target_link_libraries(lens_tests PRIVATE lenscore lensdata doctest::doctest)
target_include_directories(lens_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# Used by test_lensfile.cpp to load lensdata/cooke-s4-32mm.lens by an
# absolute path: ctest's working directory is the build tree, not the
# source tree, so a path relative to CMAKE_SOURCE_DIR would not resolve.
target_compile_definitions(lens_tests PRIVATE LENS_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
add_test(NAME lens_tests COMMAND lens_tests)
