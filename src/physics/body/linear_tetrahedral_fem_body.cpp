#include "sbs/physics/body/linear_tetrahedral_fem_body.h"

#include "sbs/physics/simulation.h"

#include <map>

namespace sbs {
namespace physics {
namespace body {

linear_tetrahedral_fem_body_t::linear_tetrahedral_fem_body_t(
    simulation_t& simulation,
    index_type const id,
    common::geometry_t const& geometry)
    : base_type(simulation, id), mechanical_model_(), visual_model_(), collision_model_()
{
    assert(geometry.geometry_type == common::geometry_t::geometry_type_t::tetrahedron);

    std::vector<Eigen::Vector3d> const tet_points = common::to_points(geometry);
    std::vector<index_type> const tet_indices     = common::to_indices(geometry);

    geometry::tetrahedral_domain_t const domain(tet_points, tet_indices);
    mechanical_model_ = mechanics::linear_tetrahedral_fem_model_t(domain);

    auto const [triangle_points, triangle_indices] = geometry::boundary_surface(domain);
    visual_model_                                  = visual::tetrahedral_fem_embedded_surface(
        triangle_points,
        triangle_indices,
        &mechanical_model_);

    collision_model_      = collision::point_bvh_model_t(&visual_model_);
    collision_model_.id() = id;
}

body_t::visual_model_type const& linear_tetrahedral_fem_body_t::visual_model() const
{
    return visual_model_;
}

body_t::collision_model_type const&
sbs::physics::body::linear_tetrahedral_fem_body_t::collision_model() const
{
    return collision_model_;
}

body_t::visual_model_type& sbs::physics::body::linear_tetrahedral_fem_body_t::visual_model()
{
    return visual_model_;
}

body_t::collision_model_type& sbs::physics::body::linear_tetrahedral_fem_body_t::collision_model()
{
    return collision_model_;
}

void sbs::physics::body::linear_tetrahedral_fem_body_t::update_visual_model()
{
    visual_model_.update();
}

void sbs::physics::body::linear_tetrahedral_fem_body_t::update_collision_model()
{
    collision_model_.update(this->simulation());
}

void sbs::physics::body::linear_tetrahedral_fem_body_t::update_physical_model()
{
    auto const& sim       = this->simulation();
    auto const idx        = this->id();
    auto const& particles = sim.particles()[idx];
    assert(particles.size() == mechanical_model_.dof_count());
    for (auto i = 0u; i < particles.size(); ++i)
    {
        mechanical_model_.dof(i) = particles[i].x();
    }
}

} // namespace body
} // namespace physics
} // namespace sbs
