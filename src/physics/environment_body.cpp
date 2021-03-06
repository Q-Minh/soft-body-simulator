#include "..\..\include\sbs\physics\environment_body.h"

#include <Discregrid/geometry/mesh_distance.hpp>
#include <Discregrid/mesh/triangle_mesh.hpp>
#include <cassert>
#include <sbs/common/geometry.h>
#include <sbs/physics/environment_body.h>

namespace sbs {
namespace physics {

environment_body_t::environment_body_t(
    simulation_t& simulation,
    index_type id,
    common::geometry_t const& geometry,
    Eigen::AlignedBox3d const& domain,
    std::array<unsigned int, 3u> const& resolution)
    : body_t(simulation, id),
      visual_model_(geometry),
      collision_model_(
          Eigen::AlignedBox3d{},
          {2, 2, 2}) /*dummy values, because Discregrid::CubicLagrangeGrid does not have a default
                        constructor*/
{
    assert(geometry.geometry_type == common::geometry_t::geometry_type_t::triangle);
    assert(geometry.has_positions());
    assert(geometry.has_indices());

    std::vector<Eigen::Vector3d> vertices{};
    std::vector<std::array<unsigned int, 3>> faces{};

    vertices.reserve(geometry.positions.size() / 3u);
    faces.reserve(geometry.indices.size() / 3u);

    for (std::size_t i = 0u; i < geometry.positions.size(); i += 3u)
    {
        auto const x = static_cast<scalar_type>(geometry.positions[i]);
        auto const y = static_cast<scalar_type>(geometry.positions[i + 1u]);
        auto const z = static_cast<scalar_type>(geometry.positions[i + 2u]);

        vertices.push_back(Eigen::Vector3d{x, y, z});
    }
    for (std::size_t i = 0u; i < geometry.indices.size(); i += 3u)
    {
        auto const v1 = static_cast<unsigned int>(geometry.indices[i]);
        auto const v2 = static_cast<unsigned int>(geometry.indices[i + 1u]);
        auto const v3 = static_cast<unsigned int>(geometry.indices[i + 2u]);

        faces.push_back({v1, v2, v3});
    }

    Discregrid::TriangleMesh mesh(vertices, faces);
    Eigen::AlignedBox3d extended_domain = domain;
    for (auto const& x : mesh.vertices())
    {
        for (auto const& x : mesh.vertices())
        {
            extended_domain.extend(x);
        }
        extended_domain.max() +=
            1.0e-3 * extended_domain.diagonal().norm() * Eigen::Vector3d::Ones();
        extended_domain.min() -=
            1.0e-3 * extended_domain.diagonal().norm() * Eigen::Vector3d::Ones();
    }

    Discregrid::MeshDistance md(mesh);

    auto const sdf = [&md](Eigen::Vector3d const& xi) {
        return md.signedDistanceCached(xi);
    };

    Discregrid::CubicLagrangeDiscreteGrid grid(extended_domain, resolution);
    grid.addFunction(sdf);

    collision_model_          = collision::sdf_model_t(grid);
    collision_model_.volume() = extended_domain;
    collision_model_.id()     = this->id();
}

environment_body_t::environment_body_t(
    simulation_t& simulation,
    index_type id,
    common::geometry_t const& geometry,
    collision::sdf_model_t const& sdf_model)
    : body_t(simulation, id), visual_model_(geometry), collision_model_(sdf_model)
{
    collision_model_.id() = this->id();
}

body_t::visual_model_type const& environment_body_t::visual_model() const
{
    return visual_model_;
}
body_t::collision_model_type const& environment_body_t::collision_model() const
{
    return collision_model_;
}
body_t::visual_model_type& environment_body_t::visual_model()
{
    return visual_model_;
}
body_t::collision_model_type& environment_body_t::collision_model()
{
    return collision_model_;
}
void environment_body_t::update_visual_model()
{
    // no-op, static mesh, but might be changed to include local transformation later
}
void environment_body_t::update_collision_model()
{
    // no-op, but might be changed to include local transformation later
}
void environment_body_t::update_physical_model()
{
    // no-op
}
void environment_body_t::transform(Eigen::Affine3d const& affine)
{
    // no-op
}

collision::sdf_model_t const& environment_body_t::sdf() const
{
    return collision_model_;
}

} // namespace physics
} // namespace sbs