add_library(lenscore INTERFACE)
target_include_directories(lenscore INTERFACE ${CMAKE_SOURCE_DIR})
target_compile_features(lenscore INTERFACE cxx_std_20)
