cmake_minimum_required(VERSION 3.14)
project(sbs VERSION 0.0.1 LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

add_executable(sbs-viewer)
set_target_properties(sbs-viewer PROPERTIES FOLDER sbs-viewer)

include(FetchContent)
FetchContent_Declare(
  glfw
  GIT_REPOSITORY https://github.com/glfw/glfw.git
  GIT_TAG        master
)
FetchContent_MakeAvailable(glfw)

add_library(glad)
target_sources(glad 
PRIVATE 
    "${CMAKE_CURRENT_LIST_DIR}/glad/src/glad.c"
)
target_include_directories(glad PUBLIC "${CMAKE_CURRENT_LIST_DIR}/glad/include")

find_package(nlohmann_json 3 CONFIG REQUIRED)
find_package(Eigen3 3.3 CONFIG REQUIRED)
find_package(glm 0.9.9 CONFIG REQUIRED)

target_compile_features(sbs-viewer PRIVATE cxx_std_17)
target_sources(sbs-viewer 
PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/main.cpp"

    # common
    "${CMAKE_CURRENT_LIST_DIR}/include/common/geometry.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/common/mesh.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/common/node.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/common/scene.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/common/mesh.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/common/node.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/common/scene.cpp"

    #geometry 
    "${CMAKE_CURRENT_LIST_DIR}/include/geometry/get_simple_bar_model.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/geometry/get_simple_cloth_model.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/geometry/get_simple_bar_model.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/geometry/get_simple_cloth_model.cpp"

    # io
    "${CMAKE_CURRENT_LIST_DIR}/include/io/endianness.hpp"
    "${CMAKE_CURRENT_LIST_DIR}/include/io/load_scene.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/io/ply.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/io/load_scene.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/io/ply.cpp"

    # physics
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/xpbd/constraint.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/xpbd/distance_constraint.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/xpbd/green_constraint.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/xpbd/solver.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/physics/xpbd/constraint.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/physics/xpbd/distance_constraint.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/physics/xpbd/green_constraint.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/physics/xpbd/solver.cpp"

    # rendering
    "${CMAKE_CURRENT_LIST_DIR}/include/rendering/camera.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/rendering/renderer.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/rendering/shader.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/rendering/camera.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/rendering/renderer.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/rendering/shader.cpp"
)
target_include_directories(sbs-viewer 
PRIVATE 
    "${CMAKE_CURRENT_LIST_DIR}/include"
)

target_link_libraries(sbs-viewer 
PRIVATE 
    glad 
    glfw 
    nlohmann_json::nlohmann_json 
    Eigen3::Eigen 
    glm::glm
)