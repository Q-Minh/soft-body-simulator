#ifndef SBS_PHYSICS_XPBD_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_CONSTRAINT_H

#include "sbs/aliases.h"

#include <Eigen/Core>
#include <vector>

namespace sbs {
namespace physics {
namespace xpbd {

// Forward declares
class simulation_t;

class constraint_t
{
  public:
    constraint_t(scalar_type alpha, scalar_type beta);
    constraint_t(
        scalar_type alpha,
        scalar_type beta,
        std::vector<index_type> const& js,
        std::vector<index_type> const& bis);

    void prepare_for_projection(simulation_t& simulation);
    virtual void project_positions(simulation_t& simulation, scalar_type dt) = 0;

    void project_positions_with_dampling(
        simulation_t& simulation,
        scalar_type const C,
        std::vector<Eigen::Vector3d> const& gradC,
        scalar_type dt);

    /**
     * @brief
     * Debug method to inspect what the values for the position corrections are
     * in subclasses.
     * @param simulation
     * @param C
     * @param gradC
     * @param dt
     * @return
     */
    std::vector<Eigen::Vector3d> get_position_corrections(
        simulation_t& simulation,
        scalar_type const C,
        std::vector<Eigen::Vector3d> const& gradC,
        scalar_type dt);

    virtual ~constraint_t() = default;

    scalar_type alpha() const;
    scalar_type beta() const;
    scalar_type lambda() const;

    scalar_type compliance() const;
    scalar_type damping() const;

    /**
     * @brief
     * Get indices of each particle involved in this constraint
     * @return
     */
    std::vector<index_type> const& js() const;
    /**
     * @brief
     * Get indices of the bodies associated with each particle involved
     * in this constraint
     * @return
     */
    std::vector<index_type> const& bs() const;

    std::vector<index_type>& js();
    std::vector<index_type>& bs();

    void set_indices(std::vector<index_type> const& js);
    void set_bodies(std::vector<index_type> const& bis);

  protected:
    virtual void prepare_for_projection_impl(simulation_t& simulation) {}

    scalar_type alpha_;           ///< Compliance coefficient
    scalar_type beta_;            ///< Damping coefficient
    scalar_type lagrange_;        ///< Lagrange multiplier
    std::vector<index_type> js_;  ///< Particle indices
    std::vector<index_type> bis_; ///< Body indices
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_CONSTRAINT_H