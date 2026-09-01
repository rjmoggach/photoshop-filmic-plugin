cmake_minimum_required(VERSION 3.24)
project(lens LANGUAGES CXX VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LENS_BUILD_TESTS "Build the test suite" ON)

add_subdirectory(lenscore)

if(LENS_BUILD_TESTS)
  include(FetchContent)
  FetchContent_Declare(doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.11)
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5)  # doctest v2.4.11 predates CMake 4's floor
  FetchContent_MakeAvailable(doctest)
  enable_testing()
  add_subdirectory(tests)
endif()
