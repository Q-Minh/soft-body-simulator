#ifndef SBS_PHYSICS_TETRAHEDRAL_BODY_H
#define SBS_PHYSICS_TETRAHEDRAL_BODY_H

#include <sbs/physics/body.h>
#include <sbs/physics/collision/bvh_model.h>
#include <sbs/physics/mesh.h>
#include <sbs/physics/tetrahedral_mesh_boundary.h>

namespace sbs {
namespace common {
class geometry_t;
} // namespace common
namespace physics {

class tetrahedral_body_t : public body_t
{
  public:
    using visual_model_type    = body_t::visual_model_type;
    using collision_model_type = body_t::collision_model_type;

    tetrahedral_body_t(std::vector<particle_t> const& particles, tetrahedron_set_t const& topology);
    tetrahedral_body_t(common::geometry_t const& geometry);

    virtual visual_model_type const& visual_model() const override;
    virtual collision_model_type const& collision_model() const override;
    virtual void update_visual_model(simulation_t const& simulation) override;
    virtual void update_collision_model(simulation_t const& simulation) override;
    virtual void update_physical_model(simulation_t const& simulation) override;

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