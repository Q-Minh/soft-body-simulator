add_executable(stvk-fem-beam)
target_sources(stvk-fem-beam PRIVATE "stvk_fem_beam.cpp")
target_link_libraries(stvk-fem-beam PRIVATE sbs::sbs)

add_executable(stvk-sph-beam)
target_sources(stvk-sph-beam PRIVATE "stvk_sph_beam.cpp")
target_link_libraries(stvk-sph-beam PRIVATE sbs::sbs)

add_executable(stvk-mls-efg-meshless-beam)
target_sources(stvk-mls-efg-meshless-beam PRIVATE "stvk_mls_efg_meshless_beam.cpp")
target_link_libraries(stvk-mls-efg-meshless-beam PRIVATE sbs::sbs)

add_executable(mls-shape-function-approximation-error)
target_sources(mls-shape-function-approximation-error PRIVATE "mls_shape_function_approximation_error.cpp")
target_link_libraries(mls-shape-function-approximation-error PRIVATE sbs::sbs)

add_executable(sandbox)
target_sources(sandbox PRIVATE "sandbox.cpp")
target_link_libraries(sandbox PRIVATE sbs::sbs)
