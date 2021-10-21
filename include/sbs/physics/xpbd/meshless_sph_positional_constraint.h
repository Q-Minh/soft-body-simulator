#ifndef SBS_PHYSICS_XPBD_MESHLESS_SPH_POSITIONAL_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_MESHLESS_SPH_POSITIONAL_CONSTRAINT_H

#include <Eigen/Core>
#include <sbs/physics/constraint.h>

namespace sbs {
namespace physics {
namespace mechanics {

class meshless_sph_body_t;
class meshless_sph_surface_vertex_t;

} // namespace mechanics
namespace xpbd {

class meshless_sph_positional_constraint_t : public constraint_t
{
  public:
    meshless_sph_positional_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        index_type bi,
        mechanics::meshless_sph_body_t const* b,
        mechanics::meshless_sph_surface_vertex_t const& vk,
        Eigen::Vector3d const& vkp,
        scalar_type mu,
        Eigen::Vector3d const& pi);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    index_type bi_;                                      ///< Body index
    mechanics::meshless_sph_body_t const* b_;            ///< The meshless bodyed
    mechanics::meshless_sph_surface_vertex_t const& vk_; ///< Surface vertex
    Eigen::Vector3d vkp_;                                ///< Surface vertex position
    scalar_type mu_;                                     ///< Penalty strength
    Eigen::Vector3d pi_;                                 ///< Target position
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_MESHLESS_SPH_POSITIONAL_CONSTRAINT_H