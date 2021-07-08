cmake_minimum_required(VERSION 3.14)
project(sbs VERSION 0.0.1 LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

if (NOT EXISTS "thirdparty/imgui")
    execute_process(
        COMMAND "git" "clone" "https://github.com/ocornut/imgui" 
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/thirdparty/" 
        RESULT_VARIABLE _imgui_clone_failure
    )
endif()

if (_imgui_clone_failure)
    message(FATAL_ERROR "Failed to download imgui")
endif()

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
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/glad/src/glad.c"
)
target_include_directories(glad PUBLIC "${CMAKE_CURRENT_LIST_DIR}/thirdparty/glad/include")

find_package(nlohmann_json 3 CONFIG REQUIRED)
find_package(Eigen3 3.3 CONFIG REQUIRED)
find_package(glm 0.9.9 CONFIG REQUIRED)

add_library(imgui)
set_target_properties(imgui PROPERTIES FOLDER imgui)
target_sources(imgui 
PRIVATE 
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/imconfig.h"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/imgui_demo.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/imgui_draw.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/imgui_internal.h"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/imgui_tables.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/imgui_widgets.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/imgui.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/imgui.h"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/imstb_rectpack.h"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/imstb_textedit.h"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/imstb_truetype.h"

    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/backends/imgui_impl_glfw.h"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/backends/imgui_impl_glfw.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/backends/imgui_impl_opengl3.h"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui/backends/imgui_impl_opengl3.cpp"
)
target_include_directories(imgui 
PUBLIC 
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty/imgui"
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

add_executable(sbs-viewer)
set_target_properties(sbs-viewer PROPERTIES FOLDER sbs-viewer)
target_compile_features(sbs-viewer PRIVATE cxx_std_17)
target_sources(sbs-viewer 
PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/main.cpp"

    # common
    "${CMAKE_CURRENT_LIST_DIR}/include/common/geometry.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/common/mesh.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/common/node.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/common/primitive.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/common/scene.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/common/geometry.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/common/mesh.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/common/node.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/common/primitive.cpp"
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
    "${CMAKE_CURRENT_LIST_DIR}/include/io/tokenize.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/io/load_scene.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/io/ply.cpp"

    # physics
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/mesh.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/physics/mesh.cpp"

    # physics/xpbd
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/xpbd/collision_constraint.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/xpbd/constraint.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/xpbd/distance_constraint.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/xpbd/green_constraint.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/xpbd/solver.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/physics/xpbd/collision_constraint.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/physics/xpbd/constraint.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/physics/xpbd/distance_constraint.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/physics/xpbd/green_constraint.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/physics/xpbd/solver.cpp"

    # physics/collision
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/collision/intersections.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/physics/collision/intersections.cpp"

    # physics/cutting
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/cutting/cut_tetrahedron.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/physics/cutting/virtual_scalpel.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/physics/cutting/cut_tetrahedron.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/physics/cutting/virtual_scalpel.cpp"

    # rendering
    "${CMAKE_CURRENT_LIST_DIR}/include/rendering/camera.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/rendering/renderer.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/rendering/shader.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/rendering/trackball_rotation_adapter.h"
    "${CMAKE_CURRENT_LIST_DIR}/include/rendering/pick.h"

    "${CMAKE_CURRENT_LIST_DIR}/src/rendering/camera.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/rendering/renderer.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/rendering/shader.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/rendering/trackball_rotation_adapter.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/rendering/pick.cpp"
)
target_include_directories(sbs-viewer 
PRIVATE 
    "${CMAKE_CURRENT_LIST_DIR}/include"
    "${CMAKE_CURRENT_LIST_DIR}/thirdparty"
)

target_link_libraries(sbs-viewer 
PRIVATE 
    glad 
    glfw 
    imgui 
    nlohmann_json::nlohmann_json 
    Eigen3::Eigen 
    glm::glm
)