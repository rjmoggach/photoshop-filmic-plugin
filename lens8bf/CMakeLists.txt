# Filmic.plugin -- a Photoshop filter (8BF), under Filter > Moggach > Filmic.
#
# The Adobe Photoshop SDK is licence-restricted and is NOT vendored in this
# repository. Point LENS_PS_SDK_ROOT at your own copy of pluginsdk/.
set(LENS_PS_SDK_ROOT "" CACHE PATH "Path to the Photoshop SDK 'pluginsdk' directory")

if(NOT LENS_PS_SDK_ROOT OR NOT EXISTS "${LENS_PS_SDK_ROOT}/photoshopapi/photoshop/PIFilter.h")
    message(STATUS "lens8bf: LENS_PS_SDK_ROOT not set to a Photoshop SDK; skipping the filter plug-in")
    return()
endif()

add_library(lens8bf MODULE
    entry.cpp
    render.cpp
    dialog.mm
)

target_compile_features(lens8bf PRIVATE cxx_std_20)

# The dialog is Objective-C++ and holds on to Cocoa objects; let ARC own their
# lifetimes rather than hand-rolling retain/release around a modal loop.
set_source_files_properties(dialog.mm PROPERTIES COMPILE_OPTIONS "-fobjc-arc")
target_link_libraries(lens8bf PRIVATE lenscore)

target_include_directories(lens8bf PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}
    "${LENS_PS_SDK_ROOT}/photoshopapi/photoshop"
    "${LENS_PS_SDK_ROOT}/photoshopapi/pica_sp"
    "${LENS_PS_SDK_ROOT}/photoshopapi/general"
    "${LENS_PS_SDK_ROOT}/samplecode/common/includes"
)

# 8BFM is what Photoshop looks for; the bundle extension must be .plugin.
set_target_properties(lens8bf PROPERTIES
    BUNDLE TRUE
    BUNDLE_EXTENSION "plugin"
    OUTPUT_NAME "Filmic"
    MACOSX_BUNDLE_INFO_PLIST "${CMAKE_CURRENT_SOURCE_DIR}/Info.plist"
    OSX_ARCHITECTURES "arm64;x86_64"
)

target_link_libraries(lens8bf PRIVATE
    "-framework Cocoa"
    "-framework CoreFoundation"
    "-framework UniformTypeIdentifiers"
)

# Photoshop reads PiPLs.json straight out of Resources -- no Rez step needed --
# and the spectral table rides along beside it.
set(_res "$<TARGET_BUNDLE_CONTENT_DIR:lens8bf>/Resources")
add_custom_command(TARGET lens8bf POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${_res}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/PiPLs.json" "${_res}/PiPLs.json"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/lensdata/rgb2spec/rec2020.bin" "${_res}/rec2020.bin"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/Filmic.icns" "${_res}/Filmic.icns"
    COMMENT "Staging PiPLs.json and the spectral table into Filmic.plugin")

# Photoshop refuses to load an unsigned or linker-signed bundle -- the same trap
# that silently stops a bundle loading. Sign it explicitly.
# Ad-hoc by default, which is all Photoshop needs: it runs with the hardened
# runtime but sets com.apple.security.cs.disable-library-validation, so an
# unsigned-by-Apple plug-in loads. Set this to a Developer ID to distribute.
set(LENS_CODESIGN_IDENTITY "-" CACHE STRING "codesign identity for the filter bundle")
add_custom_command(TARGET lens8bf POST_BUILD
    COMMAND codesign --force --sign "${LENS_CODESIGN_IDENTITY}" --timestamp=none
            "$<TARGET_BUNDLE_DIR:lens8bf>"
    COMMENT "Codesigning Filmic.plugin")

# `cmake --build build --target install-filter` drops the bundle where Photoshop
# looks for it. The app's own Plug-ins folder is root-owned, so this uses the
# per-user one, which needs no admin rights and survives Photoshop updates.
set(LENS_FILTER_INSTALL_DIR "$ENV{HOME}/Documents/Photoshop Plug-ins"
    CACHE PATH "Where to install Filmic.plugin")
add_custom_target(install-filter
    COMMAND ${CMAKE_COMMAND} -E make_directory "${LENS_FILTER_INSTALL_DIR}"
    COMMAND ${CMAKE_COMMAND} -E rm -rf "${LENS_FILTER_INSTALL_DIR}/Filmic.plugin"
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            "$<TARGET_BUNDLE_DIR:lens8bf>" "${LENS_FILTER_INSTALL_DIR}/Filmic.plugin"
    COMMAND codesign --force --sign "${LENS_CODESIGN_IDENTITY}" --timestamp=none
            "${LENS_FILTER_INSTALL_DIR}/Filmic.plugin"
    DEPENDS lens8bf
    COMMENT "Installing Filmic.plugin to ${LENS_FILTER_INSTALL_DIR}")
