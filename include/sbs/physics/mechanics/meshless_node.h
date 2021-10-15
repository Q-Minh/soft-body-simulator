#ifndef SBS_PHYSICS_MECHANICS_MESHLESS_NODE_H
#define SBS_PHYSICS_MECHANICS_MESHLESS_NODE_H

#include <Eigen/Core>
#include <sbs/aliases.h>
#include <sbs/physics/functions/kernel.h>
#include <vector>

namespace sbs {
namespace physics {

class particle_t;

namespace functions {
class kernel_t;
} // namespace functions

namespace mechanics {

class meshless_node_t
{
  public:
    meshless_node_t() = default;
    meshless_node_t(index_type const i, functions::poly6_kernel_t const& kernel);
    meshless_node_t(
        index_type const i,
        std::vector<Eigen::Vector3d const*> const& pj,
        std::vector<index_type> const& neighbours,
        functions::poly6_kernel_t const& kernel);

    index_type Ni() const;
    std::vector<index_type> const& neighbours() const;
    std::vector<scalar_type> const& Wij() const;
    std::vector<Eigen::Vector3d> const& gradWij() const;
    std::vector<scalar_type> const& Vjs() const;
    scalar_type Vi() const;
    Eigen::Matrix3d const& Li() const;
    functions::poly6_kernel_t const& kernel() const;
    Eigen::Vector3d const& Xi() const;

    void cache_Li_Vj(std::vector<meshless_node_t const*> const& neighbours);

  private:
    index_type ni_;
    std::vector<index_type> neighbours_;
    std::vector<scalar_type> Wij_;
    std::vector<Eigen::Vector3d> gradWij_;
    std::vector<scalar_type> Vjs_;
    scalar_type Vi_;
    Eigen::Matrix3d Ci_;
    functions::poly6_kernel_t kernel_;
};

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_MESHLESS_NODE_H