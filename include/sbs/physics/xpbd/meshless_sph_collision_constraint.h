#ifndef SBS_PHYSICS_XPBD_MESHLESS_SPH_COLLISION_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_MESHLESS_SPH_COLLISION_CONSTRAINT_H

#include "sbs/physics/xpbd/constraint.h"

#include <Eigen/Core>

namespace sbs {
namespace physics {
namespace mechanics {

class meshless_sph_surface_vertex_t;
class meshless_sph_body_t;

} // namespace mechanics
namespace xpbd {

class meshless_sph_collision_constraint_t : public constraint_t
{
  public:
    meshless_sph_collision_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        index_type const bi,
        mechanics::meshless_sph_body_t const* b,
        mechanics::meshless_sph_surface_vertex_t const& vk,
        Eigen::Vector3d const& vi,
        Eigen::Vector3d const& c,
        Eigen::Vector3d const& n);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    index_type bi_;                                      ///< Body index
    mechanics::meshless_sph_body_t const* b_;            ///< The meshless body
    mechanics::meshless_sph_surface_vertex_t const& vk_; ///< Penetrating vertex
    Eigen::Vector3d vi_;                                 ///< Penetrating vertex position
    Eigen::Vector3d n_;                                  ///< Contact normal
    Eigen::Vector3d c_;                                  ///< Contact point
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_MESHLESS_SPH_COLLISION_CONSTRAINT_H