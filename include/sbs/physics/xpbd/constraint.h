#ifndef SBS_PHYSICS_XPBD_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_CONSTRAINT_H

#include <memory>
#include <vector>

namespace sbs {
namespace physics {
namespace xpbd {

// forward declares
class tetrahedral_mesh_t;

class constraint_t
{
  public:
    using scalar_type = double;

    constraint_t(scalar_type const alpha) : alpha_(alpha), on_{true} {}

    virtual void project(
        std::vector<std::shared_ptr<tetrahedral_mesh_t>> const& bodies,
        scalar_type& lagrange_multiplier,
        scalar_type const dt) const = 0;

    scalar_type const& alpha() const;
    scalar_type& alpha();

    void set_as_active();
    void set_as_inactive();
    bool is_active();
    void toggle();

  private:
    scalar_type alpha_;
    bool on_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_CONSTRAINT_H