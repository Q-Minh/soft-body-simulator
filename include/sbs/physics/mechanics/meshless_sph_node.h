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

class meshless_sph_node_t
{
  public:
    meshless_sph_node_t() = default;
    meshless_sph_node_t(index_type const i, functions::poly6_kernel_t const& kernel);
    meshless_sph_node_t(
        index_type const i,
        Eigen::Vector3d const& xi,
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
    Eigen::Matrix3d const& Fi() const;
    Eigen::Matrix3d& Fi();
    Eigen::Vector3d const& xi() const;
    Eigen::Vector3d& xi();

    /**
     * The following properties represent quantities valid at time n. The sph meshless stvk
     * constraint also has duplicated logic to compute these quantities, but at every in-between
     * solve iteration. As such, the following quantities should be considered valid only either
     * before or after a solver loop.
     */

    /**
     * @brief
     * @return Green strain tensor
     */
    Eigen::Matrix3d Ei() const;
    /**
     * @brief
     * @param Ei
     * @param mu Lamé parameter
     * @param lambda Lamé parameter
     * @return The strain energy density and the stress as a pair (Psi, dPsidFi)
     */
    std::pair<scalar_type, Eigen::Matrix3d>
    Psi_dPsidFi(Eigen::Matrix3d const& Ei, scalar_type const mu, scalar_type const lambda) const;
    /**
     * @brief
     * @param Psi
     * @return Contribution of this nodal shape function to the total strain energy
     */
    scalar_type const E(scalar_type const Psi) const;
    /**
     * @brief
     * @param dPsidFi
     * @return All non-zero gradients of this nodal shape function's contribution to the total
     * strain energy w.r.t. all implicated nodes xk \in neighbours_of(ni) (this includes ni, this
     * shape function, itself).
     */
    std::vector<Eigen::Vector3d> dEdxk(Eigen::Matrix3d const& dPsidFi) const;

    void cache_Li_Vj(std::vector<meshless_sph_node_t const*> const& neighbours);

  private:
    index_type ni_;
    std::vector<index_type> neighbours_;
    std::vector<scalar_type> Wij_;
    std::vector<Eigen::Vector3d> gradWij_;
    std::vector<scalar_type> Vjs_;
    scalar_type Vi_;
    Eigen::Matrix3d Ci_;
    functions::poly6_kernel_t kernel_;
    Eigen::Matrix3d Fi_;
    Eigen::Vector3d xi_; ///< This is the interpolation coefficient, not the actual world space
                         ///< position of the meshless
};

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_MESHLESS_NODE_H