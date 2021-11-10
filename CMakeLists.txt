cmake_minimum_required(VERSION 3.14)
project(sbs VERSION 0.0.1)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Download imgui
if (NOT EXISTS "thirdparty/imgui")
    execute_process(
        COMMAND "git" "clone" "https://github.com/ocornut/imgui" 
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/" 
        RESULT_VARIABLE _imgui_clone_failure
    )
endif()

if (_imgui_clone_failure)
    message(STATUS "Failed to download imgui")
endif()

# Download implot
if (NOT EXISTS "thirdparty/implot")
    execute_process(
        COMMAND "git" "clone" "https://github.com/epezent/implot"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/" 
        RESULT_VARIABLE _implot_clone_failure
    )
endif()

if (_implot_clone_failure)
    message(STATUS "Failed to download implot")
endif()

# Download and create library target for glfw3
include(FetchContent)
FetchContent_Declare(
  glfw
  GIT_REPOSITORY https://github.com/glfw/glfw.git
  GIT_TAG        master
)
FetchContent_MakeAvailable(glfw)
set_target_properties(glfw
PROPERTIES
    DEBUG_POSTFIX "_d"
)

# Create library target for glad
add_library(glad)
set_target_properties(glad 
PROPERTIES
    DEBUG_POSTFIX "_d"
)
target_sources(glad 
PRIVATE 
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/glad/src/glad.c"
)
target_include_directories(glad 
PUBLIC 
    $<INSTALL_INTERFACE:include/glad/include>
    $<BUILD_INTERFACE:
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/glad/include
        >
)

# Get dependencies
find_package(nlohmann_json 3 CONFIG REQUIRED)
find_package(Eigen3 3.3 CONFIG REQUIRED)
find_package(glm 0.9.9 CONFIG REQUIRED)
find_package(Threads REQUIRED)

# Create library target for imgui
add_library(imgui)
set_target_properties(imgui 
PROPERTIES 
    FOLDER imgui
    DEBUG_POSTFIX "_d"    
)
target_sources(imgui 
PRIVATE 
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/imconfig.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/imgui_demo.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/imgui_draw.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/imgui_internal.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/imgui_tables.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/imgui_widgets.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/imgui.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/imgui.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/imstb_rectpack.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/imstb_textedit.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/imstb_truetype.h"

    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/backends/imgui_impl_glfw.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/backends/imgui_impl_glfw.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/backends/imgui_impl_opengl3.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/backends/imgui_impl_opengl3.cpp"
)
target_include_directories(imgui 
PUBLIC 
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui
        >
)
target_link_libraries(imgui 
PUBLIC
    glad
    glfw 
)
target_compile_definitions(imgui 
PUBLIC 
    IMGUI_IMPL_OPENGL_LOADER_GLAD=1
)

# Create library target for implot
add_library(implot)
set_target_properties(implot 
PROPERTIES 
    FOLDER implot
    DEBUG_POSTFIX "_d"    
)
target_sources(implot 
PRIVATE 
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/implot/implot.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/implot/implot.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/implot/implot_demo.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/implot/implot_internal.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/implot/implot_items.cpp"
)
target_include_directories(implot 
PUBLIC 
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/implot
        >
)
target_link_libraries(implot 
PUBLIC
    imgui
)

FetchContent_Declare(
  _discregrid
  GIT_REPOSITORY https://github.com/Q-Minh/Discregrid
  GIT_TAG        master
)
FetchContent_MakeAvailable(_discregrid)
set_target_properties(Discregrid
PROPERTIES
    DEBUG_POSTFIX "_d"
)

find_package(autodiff CONFIG REQUIRED)

# Should make the library shared in the future
add_library(sbs STATIC)
add_library(sbs::sbs ALIAS sbs)
set_target_properties(sbs 
PROPERTIES 
    FOLDER sbs
    WINDOWS_EXPORT_ALL_SYMBOLS ON 
    DEBUG_POSTFIX "_d"
)
target_compile_features(sbs PUBLIC cxx_std_17)
target_link_libraries(sbs 
PRIVATE
    nlohmann_json::nlohmann_json
PUBLIC 
    imgui 
    implot 
    Eigen3::Eigen 
    glm::glm
    Threads::Threads
    Discregrid 
    autodiff::autodiff
)

file(
    GLOB_RECURSE _sbs_sources 
    LIST_DIRECTORIES OFF 
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" 
    CONFIGURE_DEPENDS 
    "src/*.cpp")
file(
    GLOB_RECURSE _sbs_headers 
    LIST_DIRECTORIES OFF 
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" 
    CONFIGURE_DEPENDS 
    "include/*.h")

target_sources(sbs 
PRIVATE
    ${_sbs_headers}
    ${_sbs_sources}
)
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${_sbs_sources} ${_sbs_headers})
target_include_directories(sbs 
PUBLIC 
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:
        ${CMAKE_CURRENT_SOURCE_DIR}/include;
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty
        >
)

add_subdirectory(examples)

include(GNUInstallDirs)

install(
    TARGETS sbs glad imgui implot glfw Discregrid
    EXPORT sbs_targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    DIRECTORY include/ 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    DIRECTORY thirdparty/ 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    EXPORT sbs_targets
    FILE sbsTargets.cmake
    NAMESPACE sbs::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sbs
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/sbsConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/sbsConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/sbsConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sbs
)

install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/sbsConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/sbsConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sbs
)
