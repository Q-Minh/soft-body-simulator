#ifndef SBS_PHYSICS_BODY_MESHLESS_BODY_H
#define SBS_PHYSICS_BODY_MESHLESS_BODY_H

#include "sbs/aliases.h"
#include "sbs/common/geometry.h"
#include "sbs/geometry/grid.h"
#include "sbs/geometry/tetrahedral_domain.h"
#include "sbs/physics/body/body.h"
#include "sbs/physics/collision/bvh_model.h"
#include "sbs/physics/visual/meshless_embedded_surface.h"

#include <array>

namespace sbs {
namespace physics {
namespace body {

template <class MeshlessModelType>
class meshless_body_t : public body_t
{
  public:
    using base_type                    = body_t;
    using meshless_model_type          = MeshlessModelType;
    using self_type                    = meshless_body_t<meshless_model_type>;
    using embedded_visual_surface_type = visual::meshless_embedded_surface_t<meshless_model_type>;

    meshless_body_t(
        xpbd::simulation_t& simulation,
        index_type const id,
        common::geometry_t const& geometry,
        std::array<unsigned int, 3u> const& resolution,
        scalar_type support);

    virtual visual_model_type const& visual_model() const override;
    virtual collision_model_type const& collision_model() const override;
    virtual visual_model_type& visual_model() override;
    virtual collision_model_type& collision_model() override;
    virtual void update_visual_model() override;
    virtual void update_collision_model() override;
    virtual void update_physical_model() override;

    meshless_model_type const& get_mechanical_model() const;
    embedded_visual_surface_type const& get_visual_model() const;
    collision::point_bvh_model_t const& get_collision_model() const;

    meshless_model_type& get_mechanical_model();
    embedded_visual_surface_type& get_visual_model();
    collision::point_bvh_model_t& get_collision_model();

  private:
    meshless_model_type mechanical_model_;
    embedded_visual_surface_type visual_model_;
    collision::point_bvh_model_t collision_model_;
};

template <class MeshlessModelType>
inline meshless_body_t<MeshlessModelType>::meshless_body_t(
    xpbd::simulation_t& simulation,
    index_type const id,
    common::geometry_t const& geometry,
    std::array<unsigned int, 3u> const& resolution,
    scalar_type support)
    : base_type(simulation, id), mechanical_model_(), visual_model_(), collision_model_()
{
    assert(geometry.geometry_type == common::geometry_t::geometry_type_t::tetrahedron);

    std::vector<Eigen::Vector3d> const tet_points = common::to_points(geometry);
    std::vector<index_type> const tet_indices     = common::to_indices(geometry);

    geometry::tetrahedral_domain_t const domain(tet_points, tet_indices);
    geometry::grid_t const grid(
        domain,
        Eigen::Vector3i{
            static_cast<int>(resolution[0]),
            static_cast<int>(resolution[1]),
            static_cast<int>(resolution[2])});

    mechanical_model_ = meshless_model_type(domain, grid, support);

    auto const [triangle_points, triangle_indices] = geometry::boundary_surface(domain);
    visual_model_ =
        embedded_visual_surface_type(triangle_points, triangle_indices, &mechanical_model_);

    collision_model_      = collision::point_bvh_model_t(&visual_model_);
    collision_model_.id() = id;
}

template <class MeshlessModelType>
inline body_t::visual_model_type const& meshless_body_t<MeshlessModelType>::visual_model() const
{
    return visual_model_;
}

template <class MeshlessModelType>
inline body_t::collision_model_type const&
meshless_body_t<MeshlessModelType>::collision_model() const
{
    return collision_model_;
}

template <class MeshlessModelType>
inline body_t::visual_model_type& meshless_body_t<MeshlessModelType>::visual_model()
{
    return visual_model_;
}

template <class MeshlessModelType>
inline body_t::collision_model_type& meshless_body_t<MeshlessModelType>::collision_model()
{
    return collision_model_;
}

template <class MeshlessModelType>
inline void meshless_body_t<MeshlessModelType>::update_visual_model()
{
    visual_model_.update();
}

template <class MeshlessModelType>
inline void meshless_body_t<MeshlessModelType>::update_collision_model()
{
    collision_model_.update(this->simulation());
}

template <class MeshlessModelType>
inline void meshless_body_t<MeshlessModelType>::update_physical_model()
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

template <class MeshlessModelType>
inline MeshlessModelType const& meshless_body_t<MeshlessModelType>::get_mechanical_model() const
{
    return mechanical_model_;
}

template <class MeshlessModelType>
inline typename meshless_body_t<MeshlessModelType>::embedded_visual_surface_type const&
meshless_body_t<MeshlessModelType>::get_visual_model() const
{
    return visual_model_;
}

template <class MeshlessModelType>
inline collision::point_bvh_model_t const&
meshless_body_t<MeshlessModelType>::get_collision_model() const
{
    return collision_model_;
}

template <class MeshlessModelType>
inline MeshlessModelType& meshless_body_t<MeshlessModelType>::get_mechanical_model()
{
    return mechanical_model_;
}

template <class MeshlessModelType>
inline typename meshless_body_t<MeshlessModelType>::embedded_visual_surface_type&
meshless_body_t<MeshlessModelType>::get_visual_model()
{
    return visual_model_;
}

template <class MeshlessModelType>
inline collision::point_bvh_model_t& meshless_body_t<MeshlessModelType>::get_collision_model()
{
    return collision_model_;
}

} // namespace body
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_BODY_MESHLESS_BODY_H