add_executable(lens_tests test_smoke.cpp test_image.cpp test_transfer.cpp)
target_link_libraries(lens_tests PRIVATE lenscore doctest::doctest)
add_test(NAME lens_tests COMMAND lens_tests)
