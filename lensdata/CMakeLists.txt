include(FetchContent)
FetchContent_Declare(nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.11.3)
FetchContent_MakeAvailable(nlohmann_json)

add_library(lensdata INTERFACE)
target_include_directories(lensdata INTERFACE ${CMAKE_SOURCE_DIR})
target_link_libraries(lensdata INTERFACE lenscore nlohmann_json::nlohmann_json)
