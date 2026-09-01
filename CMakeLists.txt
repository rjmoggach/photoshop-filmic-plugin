cmake_minimum_required(VERSION 3.24)
project(lens LANGUAGES CXX VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Build type; empty default would mean no optimisation" FORCE)
endif()

option(LENS_BUILD_TESTS "Build the test suite" ON)

add_subdirectory(lenscore)
add_subdirectory(rgb2spec)

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
