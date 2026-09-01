add_executable(lens_tests test_smoke.cpp test_image.cpp test_transfer.cpp test_cie.cpp test_upsample.cpp test_spectable.cpp test_dispersion.cpp test_warp.cpp)
target_link_libraries(lens_tests PRIVATE lenscore doctest::doctest)
add_test(NAME lens_tests COMMAND lens_tests)
