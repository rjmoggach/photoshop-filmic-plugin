# The .uxpaddon is a plain shared library with its extension renamed. That is
# the whole trick, and it is why one target serves both platforms -- the SDK's
# Xcode and Visual Studio projects are samples, not a requirement.
#
# Adobe's UXP Hybrid SDK is licence-restricted and is NOT vendored here. Point
# UXP_ADDON_SDK_ROOT at your unpacked copy:
#   cmake -B build -S . -DUXP_ADDON_SDK_ROOT=/path/to/uxp-hybrid-plugin-sdk-main

if(NOT UXP_ADDON_SDK_ROOT)
  message(STATUS "lensaddon: UXP_ADDON_SDK_ROOT not set, skipping the plugin addon.")
  message(STATUS "           Set it to build the loadable Photoshop plugin.")
  return()
endif()

set(_sdk_src "${UXP_ADDON_SDK_ROOT}/src")
if(NOT EXISTS "${_sdk_src}/utilities/UxpAddon.h")
  message(FATAL_ERROR "lensaddon: no UxpAddon.h under ${_sdk_src}/utilities -- is UXP_ADDON_SDK_ROOT correct?")
endif()

add_library(lensaddon SHARED
  module.cpp
  "${_sdk_src}/utilities/UxpAddon.cpp"
  "${_sdk_src}/utilities/UxpTask.cpp"
  "${_sdk_src}/utilities/UxpValue.cpp")

target_include_directories(lensaddon PRIVATE
  "${_sdk_src}/utilities"
  "${_sdk_src}/api")

target_link_libraries(lensaddon PRIVATE lenscore)

# Photoshop loads the file by extension, not by platform convention.
set_target_properties(lensaddon PROPERTIES
  PREFIX ""
  OUTPUT_NAME "lens"
  SUFFIX ".uxpaddon"
  CXX_VISIBILITY_PRESET default)

if(APPLE)
  # Universal so one build serves Apple Silicon and Intel.
  set_target_properties(lensaddon PROPERTIES OSX_ARCHITECTURES "arm64;x86_64")
endif()

# Stage into the layout UXP resolves: <platform>/<arch>/<name>.uxpaddon.
# Our macOS build is universal, so the same file serves both arch folders --
# Adobe's own template ships separate arm64 and x64 binaries, but one fat
# binary in both places is equivalent and halves what we have to build.
if(APPLE)
  set(_stages "${CMAKE_SOURCE_DIR}/plugin/mac/arm64" "${CMAKE_SOURCE_DIR}/plugin/mac/x64")
else()
  set(_stages "${CMAKE_SOURCE_DIR}/plugin/win/x64")
endif()

add_custom_command(TARGET lensaddon POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/plugin/data"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${CMAKE_SOURCE_DIR}/lensdata/rgb2spec/rec2020.bin"
          "${CMAKE_SOURCE_DIR}/plugin/data/"
  COMMENT "Staging lens.uxpaddon and the spectral table into plugin/")

foreach(_stage IN LISTS _stages)
  add_custom_command(TARGET lensaddon POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${_stage}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:lensaddon>" "${_stage}/")
endforeach()
