add_executable(fem-beam)
target_sources(fem-beam PRIVATE "fem_beam.cpp")
target_link_libraries(fem-beam PRIVATE sbs::sbs)

add_executable(fem-beam-autodiff)
target_sources(fem-beam-autodiff PRIVATE "fem_beam_autodiff.cpp")
target_link_libraries(fem-beam-autodiff PRIVATE sbs::sbs)

add_executable(sph-beam)
target_sources(sph-beam PRIVATE "sph_beam.cpp")
target_link_libraries(sph-beam PRIVATE sbs::sbs)

add_executable(sph-efg-beam)
target_sources(sph-efg-beam PRIVATE "sph_efg_beam.cpp")
target_link_libraries(sph-efg-beam PRIVATE sbs::sbs)

add_executable(mls-efg-meshless-beam)
target_sources(mls-efg-meshless-beam PRIVATE "mls_efg_meshless_beam.cpp")
target_link_libraries(mls-efg-meshless-beam PRIVATE sbs::sbs)

add_executable(mls-shape-function-approximation-error)
target_sources(mls-shape-function-approximation-error PRIVATE "mls_shape_function_approximation_error.cpp")
target_link_libraries(mls-shape-function-approximation-error PRIVATE sbs::sbs)

add_executable(fem-sph-beam)
target_sources(fem-sph-beam PRIVATE "fem_sph_beam.cpp")
target_link_libraries(fem-sph-beam PRIVATE sbs::sbs)

add_executable(bezier-curve-swept-surface)
target_sources(bezier-curve-swept-surface PRIVATE "bezier_curve_swept_surface.cpp")
target_link_libraries(bezier-curve-swept-surface PRIVATE sbs::sbs)

add_executable(sandbox)
target_sources(sandbox PRIVATE "sandbox.cpp")
target_link_libraries(sandbox PRIVATE sbs::sbs)
