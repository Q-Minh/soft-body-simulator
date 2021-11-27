#ifndef SBS_PHYSICS_BODY_FEM_MIXED_BODY_H
#define SBS_PHYSICS_BODY_FEM_MIXED_BODY_H

#include "body.h"
#include "sbs/common/geometry.h"
#include "sbs/geometry/grid.h"
#include "sbs/geometry/tetrahedral_domain.h"
#include "sbs/physics/collision/bvh_model.h"
#include "sbs/physics/visual/fem_mixed_embedded_surface.h"

#include <Eigen/Core>

namespace sbs {
namespace physics {
namespace body {

template <class FemMixedModelType>
class fem_mixed_body_t : public body_t
{
  public:
    using base_type                    = body_t;
    using fem_mixed_model_type         = FemMixedModelType;
    using self_type                    = fem_mixed_body_t<fem_mixed_model_type>;
    using embedded_visual_surface_type = visual::fem_mixed_embedded_surface_t<fem_mixed_model_type>;

    fem_mixed_body_t(
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

    fem_mixed_model_type const& get_mechanical_model() const;
    embedded_visual_surface_type const& get_visual_model() const;
    collision::point_bvh_model_t const& get_collision_model() const;

    fem_mixed_model_type& get_mechanical_model();
    embedded_visual_surface_type& get_visual_model();
    collision::point_bvh_model_t& get_collision_model();

    void set_fem_particle_offset(std::size_t offset);
    void set_meshless_particle_offset(std::size_t offset);

  private:
    fem_mixed_model_type mechanical_model_;
    embedded_visual_surface_type visual_model_;
    collision::point_bvh_model_t collision_model_;

    std::size_t fem_particle_offset_;
    std::size_t meshless_particle_offset_;
};

template <class FemMixedModelType>
inline fem_mixed_body_t<FemMixedModelType>::fem_mixed_body_t(
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

    mechanical_model_ = fem_mixed_model_type(domain, grid, support);

    auto const [triangle_points, triangle_indices] = geometry::boundary_surface(domain);
    visual_model_ =
        embedded_visual_surface_type(triangle_points, triangle_indices, &mechanical_model_);

    collision_model_      = collision::point_bvh_model_t(&visual_model_);
    collision_model_.id() = id;
}

template <class FemMixedModelType>
inline typename body_t::visual_model_type const&
fem_mixed_body_t<FemMixedModelType>::visual_model() const
{
    return visual_model_;
}

template <class FemMixedModelType>
inline typename body_t::collision_model_type const&
fem_mixed_body_t<FemMixedModelType>::collision_model() const
{
    return collision_model_;
}

template <class FemMixedModelType>
inline typename body_t::visual_model_type& fem_mixed_body_t<FemMixedModelType>::visual_model()
{
    return visual_model_;
}

template <class FemMixedModelType>
inline typename body_t::collision_model_type& fem_mixed_body_t<FemMixedModelType>::collision_model()
{
    return collision_model_;
}

template <class FemMixedModelType>
inline void fem_mixed_body_t<FemMixedModelType>::update_visual_model()
{
    visual_model_.update();
}

template <class FemMixedModelType>
inline void fem_mixed_body_t<FemMixedModelType>::update_collision_model()
{
    collision_model_.update(this->simulation());
}

template <class FemMixedModelType>
inline void fem_mixed_body_t<FemMixedModelType>::update_physical_model()
{
    auto const& sim                   = this->simulation();
    auto const idx                    = this->id();
    auto const& particles             = sim.particles()[idx];
    auto& meshless_model              = mechanical_model_.meshless_model();
    auto const num_fem_particles      = mechanical_model_.dof_count();
    auto const num_meshless_particles = meshless_model.dof_count();
    assert(particles.size() == num_fem_particles + num_meshless_particles);

    for (auto i = 0u; i < num_fem_particles; ++i)
    {
        auto const idx           = fem_particle_offset_ + i;
        mechanical_model_.dof(i) = particles[idx].x();
    }
    for (auto j = 0u; j < num_meshless_particles; ++j)
    {
        auto const idx        = meshless_particle_offset_ + j;
        meshless_model.dof(j) = particles[idx].x();
    }
}

template <class FemMixedModelType>
inline FemMixedModelType const& fem_mixed_body_t<FemMixedModelType>::get_mechanical_model() const
{
    return mechanical_model_;
}

template <class FemMixedModelType>
inline typename fem_mixed_body_t<FemMixedModelType>::embedded_visual_surface_type const&
fem_mixed_body_t<FemMixedModelType>::get_visual_model() const
{
    return visual_model_;
}

template <class FemMixedModelType>
inline collision::point_bvh_model_t const&
fem_mixed_body_t<FemMixedModelType>::get_collision_model() const
{
    return collision_model_;
}

template <class FemMixedModelType>
inline FemMixedModelType& fem_mixed_body_t<FemMixedModelType>::get_mechanical_model()
{
    return mechanical_model_;
}

template <class FemMixedModelType>
inline typename fem_mixed_body_t<FemMixedModelType>::embedded_visual_surface_type&
fem_mixed_body_t<FemMixedModelType>::get_visual_model()
{
    return visual_model_;
}

template <class FemMixedModelType>
inline collision::point_bvh_model_t& fem_mixed_body_t<FemMixedModelType>::get_collision_model()
{
    return collision_model_;
}

template <class FemMixedModelType>
inline void fem_mixed_body_t<FemMixedModelType>::set_fem_particle_offset(std::size_t offset)
{
    fem_particle_offset_ = offset;
}

template <class FemMixedModelType>
inline void fem_mixed_body_t<FemMixedModelType>::set_meshless_particle_offset(std::size_t offset)
{
    meshless_particle_offset_ = offset;
}

} // namespace body
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_BODY_FEM_MIXED_BODY_H
