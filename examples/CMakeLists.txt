add_executable(stvk-fem-beam)
target_sources(stvk-fem-beam PRIVATE "stvk_fem_beam.cpp")
target_link_libraries(stvk-fem-beam PRIVATE sbs::sbs)

add_executable(stvk-meshless-beam)
target_sources(stvk-meshless-beam PRIVATE "stvk_meshless_beam.cpp")
target_link_libraries(stvk-meshless-beam PRIVATE sbs::sbs)

add_executable(stvk-hybrid-mesh-meshless-beam)
target_sources(stvk-hybrid-mesh-meshless-beam PRIVATE "stvk_hybrid_mesh_meshless_beam.cpp")
target_link_libraries(stvk-hybrid-mesh-meshless-beam PRIVATE sbs::sbs)

add_executable(stvk-mls-efg-meshless-beam)
target_sources(stvk-mls-efg-meshless-beam PRIVATE "stvk_mls_efg_meshless_beam.cpp")
target_link_libraries(stvk-mls-efg-meshless-beam PRIVATE sbs::sbs)

add_executable(sandbox)
target_sources(sandbox PRIVATE "sandbox.cpp")
target_link_libraries(sandbox PRIVATE sbs::sbs)
