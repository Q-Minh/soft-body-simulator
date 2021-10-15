add_executable(stvk-fem-beam)
target_sources(stvk-fem-beam PRIVATE "stvk_fem_beam.cpp")
target_link_libraries(stvk-fem-beam PRIVATE sbs::sbs)

add_executable(stvk-meshless-beam)
target_sources(stvk-meshless-beam PRIVATE "stvk_meshless_beam.cpp")
target_link_libraries(stvk-meshless-beam PRIVATE sbs::sbs)
