#ifndef SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_MLS_NODE_H
#define SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_MLS_NODE_H

#include <Eigen/Core>
#include <array>
#include <optional>
#include <sbs/physics/functions/kernel.h>

namespace sbs {
namespace physics {
namespace mechanics {

class hybrid_mesh_meshless_mls_body_t;

class hybrid_mesh_meshless_mls_node_t
{
  public:
    hybrid_mesh_meshless_mls_node_t(
        index_type const i,
        functions::poly6_kernel_t const& kernel,
        hybrid_mesh_meshless_mls_body_t& body);

    void initialize(
        Eigen::Vector3d const& xi,
        std::vector<Eigen::Vector3d const*> const& pj,
        std::vector<index_type> const& neighbours,
        index_type const ti);

    bool is_mixed_particle() const;

    index_type Ni() const;
    functions::poly6_kernel_t const& kernel() const;
    functions::poly6_kernel_t& kernel();
    hybrid_mesh_meshless_mls_body_t const& body() const;
    std::vector<index_type> const& neighbours() const;
    std::vector<scalar_type> const& Wijs() const;
    std::vector<Eigen::Vector3d> const& gradWijs() const;
    scalar_type Vi() const;
    Eigen::Vector4d const& alphaI() const;
    std::vector<scalar_type> const& phi_js() const;
    std::vector<Eigen::Vector3d> const& grad_phi_js() const;
    std::array<std::optional<scalar_type>, 4u> const& mesh_phi_js() const;
    std::array<std::optional<Eigen::Vector3d>, 4u> const& mesh_grad_phi_js() const;
    Eigen::Vector3d const& Xi() const;
    Eigen::Vector3d& Xi();
    Eigen::Matrix3d const& Fi() const;
    Eigen::Matrix3d& Fi();
    Eigen::Vector3d const& xi() const;
    Eigen::Vector3d& xi();
    Eigen::Matrix3d const& Ei() const;
    Eigen::Matrix3d& Ei();

    index_type ti() const;
    std::array<index_type, 4u> const& vis() const;

    Eigen::Vector4d polynomial(Eigen::Vector3d const& X) const;

  private:
    index_type ni_;
    functions::poly6_kernel_t kernel_;
    hybrid_mesh_meshless_mls_body_t& body_;
    std::vector<index_type> neighbours_;
    std::vector<scalar_type> Wijs_;
    std::vector<Eigen::Vector3d> gradWijs_;
    scalar_type Vi_;
    Eigen::Vector4d alpha_i_;
    Eigen::Matrix4d M_;
    Eigen::Matrix4d Minv_;
    std::array<Eigen::Matrix4d, 3u> gradM_;
    std::array<Eigen::Matrix4d, 3u> gradMinv_;
    std::vector<scalar_type> phi_js_;
    std::vector<Eigen::Vector3d> grad_phi_js_;
    index_type ti_; ///< Boundary tetrahedron in which this meshless node lives
    std::array<std::optional<scalar_type>, 4u> mesh_phi_js_;
    std::array<std::optional<Eigen::Vector3d>, 4u> mesh_grad_phi_js_;

    Eigen::Matrix3d Fi_;
    Eigen::Vector3d xi_;
    Eigen::Matrix3d Ei_;
};

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_MLS_NODE_H