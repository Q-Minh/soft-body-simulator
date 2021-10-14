#ifndef SBS_PHYSICS_TETRAHEDRAL_BODY_H
#define SBS_PHYSICS_TETRAHEDRAL_BODY_H

#include <sbs/aliases.h>
#include <sbs/physics/body.h>
#include <sbs/physics/collision/bvh_model.h>
#include <sbs/physics/particle.h>
#include <sbs/physics/tetrahedral_mesh_boundary.h>
#include <sbs/physics/topology.h>

namespace sbs {
namespace common {

struct geometry_t;

} // namespace common

namespace physics {

class simulation_t;

class tetrahedral_body_t : public body_t
{
  public:
    using visual_model_type    = body_t::visual_model_type;
    using collision_model_type = body_t::collision_model_type;

    tetrahedral_body_t(
        simulation_t& simulation,
        index_type id,
        std::vector<Eigen::Vector3d> const& positions,
        tetrahedron_set_t const& topology);
    tetrahedral_body_t(simulation_t& simulation, index_type id, common::geometry_t const& geometry);

    virtual visual_model_type const& visual_model() const override;
    virtual collision_model_type const& collision_model() const override;
    virtual visual_model_type& visual_model() override;
    virtual collision_model_type& collision_model() override;
    virtual void update_visual_model() override;
    virtual void update_collision_model() override;
    virtual void update_physical_model() override;
    virtual void transform(Eigen::Affine3d const& affine) override;

    tetrahedron_set_t const& physical_model() const;

    tetrahedral_mesh_boundary_t const& surface_mesh() const;
    tetrahedral_mesh_boundary_t& surface_mesh();

    collision::point_bvh_model_t const& bvh() const;
    collision::point_bvh_model_t& bvh();

  protected:
    void update_visual_model(std::vector<particle_t> const& particles);

  private:
    tetrahedron_set_t physical_model_;
    tetrahedral_mesh_boundary_t visual_model_;
    collision::point_bvh_model_t collision_model_;
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_TETRAHEDRAL_BODY_H