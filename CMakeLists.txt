cmake_minimum_required(VERSION 3.14)
project(sbs VERSION 0.0.1)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

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

find_package(nlohmann_json 3 CONFIG REQUIRED)
find_package(Eigen3 3.3 CONFIG REQUIRED)
find_package(glm 0.9.9 CONFIG REQUIRED)
find_package(Threads REQUIRED)

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
    Eigen3::Eigen 
    glm::glm
    Threads::Threads
    Discregrid
)
target_sources(sbs 
PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/aliases.h"

    # common
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/common/geometry.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/common/mesh.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/common/node.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/common/primitive.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/common/scene.h"

    "${CMAKE_CURRENT_SOURCE_DIR}/src/common/geometry.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/common/mesh.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/common/node.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/common/primitive.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/common/scene.cpp"

    #geometry 
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/geometry/get_simple_bar_model.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/geometry/get_simple_cloth_model.h"

    "${CMAKE_CURRENT_SOURCE_DIR}/src/geometry/get_simple_bar_model.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/geometry/get_simple_cloth_model.cpp"

    # io
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/io/endianness.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/io/load_scene.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/io/ply.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/io/tokenize.h"

    "${CMAKE_CURRENT_SOURCE_DIR}/src/io/load_scene.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/io/ply.cpp"

    # physics
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/body.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/constraint.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/environment_body.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/gauss_seidel_solver.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/particle.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/simulation.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/solver.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/tetrahedral_body.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/timestep.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/tetrahedral_mesh_boundary.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/topology.h"

    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/body.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/constraint.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/environment_body.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/gauss_seidel_solver.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/particle.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/simulation.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/solver.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/tetrahedral_body.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/timestep.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/tetrahedral_mesh_boundary.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/topology.cpp"

    # physics/xpbd
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/xpbd/collision_constraint.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/xpbd/contact_handler.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/xpbd/distance_constraint.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/xpbd/green_constraint.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/xpbd/simulation_parameters.h"

    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/xpbd/collision_constraint.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/xpbd/contact_handler.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/xpbd/distance_constraint.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/xpbd/green_constraint.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/xpbd/simulation_parameters.cpp"

    # physics/collision
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/collision/brute_force_cd_system.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/collision/brute_force_collision_detector.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/collision/bvh.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/collision/bvh_model.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/collision/cd_system.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/collision/collision_detector.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/collision/collision_model.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/collision/contact.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/collision/intersections.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/collision/sdf_model.h"

    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/collision/brute_force_cd_system.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/collision/brute_force_collision_detector.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/collision/bvh_model.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/collision/cd_system.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/collision/collision_detector.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/collision/collision_model.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/collision/contact.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/collision/intersections.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/collision/sdf_model.cpp"

    # physics/cutting
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/cutting/cut_tetrahedron.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/physics/cutting/virtual_scalpel.h"

    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/cutting/cut_tetrahedron.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/physics/cutting/virtual_scalpel.cpp"

    # rendering
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/rendering/camera.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/rendering/renderer.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/rendering/shader.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/rendering/trackball_rotation_adapter.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/rendering/physics_timestep_throttler.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/sbs/rendering/pick.h"

    "${CMAKE_CURRENT_SOURCE_DIR}/src/rendering/camera.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/rendering/renderer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/rendering/shader.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/rendering/trackball_rotation_adapter.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/rendering/physics_timestep_throttler.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/rendering/pick.cpp"
)
target_include_directories(sbs 
PUBLIC 
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:
        ${CMAKE_CURRENT_SOURCE_DIR}/include;
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty
        >
)

add_executable(tester)
target_sources(tester PRIVATE main.cpp)
target_link_libraries(tester PRIVATE sbs)

include(GNUInstallDirs)

install(
    TARGETS sbs glad imgui glfw Discregrid
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
