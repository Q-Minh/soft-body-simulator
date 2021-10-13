#include <cassert>
#include <sbs/common/geometry.h>
#include <sbs/physics/particle.h>
#include <sbs/physics/simulation.h>
#include <sbs/physics/tetrahedral_body.h>

namespace sbs {
namespace physics {

tetrahedral_body_t::tetrahedral_body_t(
    std::vector<particle_t> const& particles,
    tetrahedron_set_t const& topology)
    : physical_model_(topology), visual_model_(&physical_model_), collision_model_()
{
    update_visual_model(particles);
}

tetrahedral_body_t::tetrahedral_body_t(common::geometry_t const& geometry)
    : physical_model_(), visual_model_(), collision_model_()
{
    assert(geometry.geometry_type == common::geometry_t::geometry_type_t::tetrahedron);
    assert(geometry.has_indices());
    assert(geometry.has_positions());

    physical_model_.reserve_vertices(geometry.positions.size() / 3u);

    for (std::size_t i = 0u; i < geometry.indices.size(); i += 4u)
    {
        sbs::physics::tetrahedron_t const tetrahedron{
            static_cast<index_type>(geometry.indices[i]),
            static_cast<index_type>(geometry.indices[i + 1u]),
            static_cast<index_type>(geometry.indices[i + 2u]),
            static_cast<index_type>(geometry.indices[i + 3u])};

        physical_model_.add_tetrahedron(tetrahedron);
    }
    visual_model_ = tetrahedral_mesh_boundary_t(&physical_model_);

    for (std::size_t i = 0u; i < visual_model_.vertex_count(); ++i)
    {
        auto const particle_index = visual_model_.from_surface_vertex(i);

        auto const idx      = i * 3u;
        scalar_type const x = static_cast<scalar_type>(geometry.positions[idx]);
        scalar_type const y = static_cast<scalar_type>(geometry.positions[idx + 1u]);
        scalar_type const z = static_cast<scalar_type>(geometry.positions[idx + 2u]);

        float const r = static_cast<float>(geometry.colors[idx] / 255.f);
        float const g = static_cast<float>(geometry.colors[idx + 1u] / 255.f);
        float const b = static_cast<float>(geometry.colors[idx + 2u] / 255.f);

        visual_model_.vertex(i).position = Eigen::Vector3d{x, y, z};
        visual_model_.vertex(i).color    = Eigen::Vector3f{r, g, b};
    }
    visual_model_.compute_normals();

    collision_model_ = collision::point_bvh_model_t(&visual_model_);
}

tetrahedral_body_t::visual_model_type const& tetrahedral_body_t::visual_model() const
{
    return visual_model_;
}

tetrahedral_body_t::collision_model_type const& tetrahedral_body_t::collision_model() const
{
    return collision_model_;
}

body_t::visual_model_type& tetrahedral_body_t::visual_model()
{
    return visual_model_;
}

body_t::collision_model_type& tetrahedral_body_t::collision_model()
{
    return collision_model_;
}

void tetrahedral_body_t::update_visual_model(simulation_t const& simulation)
{
    auto const& particles = simulation.particles()[id()];
    update_visual_model(particles);
}

void tetrahedral_body_t::update_collision_model(simulation_t const& simulation)
{
    collision_model_.update(simulation);
}

void tetrahedral_body_t::update_physical_model(simulation_t const& simulation)
{
    // no-op
}

tetrahedron_set_t const& tetrahedral_body_t::physical_model() const
{
    return physical_model_;
}

tetrahedral_mesh_boundary_t const& tetrahedral_body_t::surface_mesh() const
{
    return visual_model_;
}
tetrahedral_mesh_boundary_t& tetrahedral_body_t::surface_mesh()
{
    return visual_model_;
}

collision::point_bvh_model_t const& tetrahedral_body_t::bvh() const
{
    return collision_model_;
}
collision::point_bvh_model_t& tetrahedral_body_t::bvh()
{
    return collision_model_;
}

void tetrahedral_body_t::update_visual_model(std::vector<particle_t> const& particles)
{
    for (std::size_t i = 0u; i < visual_model_.vertex_count(); ++i)
    {
        auto const particle_index        = visual_model_.from_surface_vertex(i);
        visual_model_.vertex(i).position = particles[particle_index].x();
    }
    visual_model_.compute_normals();
}

} // namespace physics
} // namespace sbs