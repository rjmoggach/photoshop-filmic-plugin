add_executable(lens_tests test_smoke.cpp test_image.cpp test_transfer.cpp test_cie.cpp test_upsample.cpp test_spectable.cpp test_dispersion.cpp test_warp.cpp test_vignette.cpp test_zernike.cpp test_pupil.cpp test_fft.cpp test_psf.cpp test_eff.cpp test_psfrings.cpp test_targets.cpp test_metrics.cpp test_pipeline.cpp test_acceptance.cpp)
target_link_libraries(lens_tests PRIVATE lenscore doctest::doctest)
target_include_directories(lens_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME lens_tests COMMAND lens_tests)
