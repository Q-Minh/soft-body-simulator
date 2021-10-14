#ifndef SBS_PHYSICS_BODY_H
#define SBS_PHYSICS_BODY_H

#include <Eigen/Geometry>
#include <memory>
#include <sbs/aliases.h>
#include <sbs/common/node.h>
#include <sbs/physics/collision/collision_model.h>

namespace sbs {
namespace physics {

// Forward declares
class simulation_t;

class body_t
{
  public:
    using visual_model_type    = common::renderable_node_t;
    using collision_model_type = collision::collision_model_t;

    body_t(simulation_t& simulation, index_type id) : simulation_(simulation), id_(id) {}

    virtual visual_model_type const& visual_model() const       = 0;
    virtual collision_model_type const& collision_model() const = 0;
    virtual visual_model_type& visual_model()                   = 0;
    virtual collision_model_type& collision_model()             = 0;
    virtual void update_visual_model()                          = 0;
    virtual void update_collision_model()                       = 0;
    virtual void update_physical_model()                        = 0;
    virtual void transform(Eigen::Affine3d const& affine)       = 0;

    index_type id() const;
    index_type& id();

    virtual ~body_t() = default;

    simulation_t const& simulation() const { return simulation_; }

  protected:
    simulation_t& simulation() { return simulation_; }

  private:
    index_type id_;
    simulation_t& simulation_;
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_BODY_H