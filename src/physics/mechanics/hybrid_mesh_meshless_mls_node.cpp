#include <sbs/physics/mechanics/hybrid_mesh_meshless_mls_body.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_mls_node.h>

namespace sbs {
namespace physics {
namespace mechanics {

hybrid_mesh_meshless_mls_node_t::hybrid_mesh_meshless_mls_node_t(
    index_type const i,
    functions::poly6_kernel_t const& kernel,
    hybrid_mesh_meshless_mls_body_t& body)
    : ni_(i),
      kernel_(kernel),
      body_(body),
      neighbours_(),
      Wijs_(),
      gradWijs_(),
      Vi_(),
      alpha_i_(),
      M_(),
      Minv_(),
      gradM_(),
      gradMinv_(),
      phi_js_(),
      grad_phi_js_(),
      ti_(std::numeric_limits<index_type>::max()),
      mesh_phi_js_(),
      mesh_grad_phi_js_(),
      Fi_(),
      xi_()
{
}

void hybrid_mesh_meshless_mls_node_t::initialize(
    Eigen::Vector3d const& xi,
    std::vector<Eigen::Vector3d const*> const& pj,
    std::vector<index_type> const& neighbours,
    index_type const ti)
{
    neighbours_ = neighbours;
    ti_         = ti;
    Fi_.setIdentity();
    xi_ = xi;
    Ei_.setZero();

    auto const& W = kernel_;

    M_.setZero();
    Minv_.setZero();
    gradM_[0].setZero();
    gradM_[1].setZero();
    gradM_[2].setZero();
    gradMinv_[0].setZero();
    gradMinv_[1].setZero();
    gradMinv_[2].setZero();
    alpha_i_.setZero();

    Vi_ = 0.;
    std::vector<Eigen::Vector4d> Ajs{};
    std::vector<std::array<Eigen::Vector4d, 3u>> gradAjs{};
    Ajs.reserve(neighbours.size());
    gradAjs.reserve(neighbours.size());
    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        index_type const j            = neighbours[a];
        auto const& neighbour         = body_.meshless_nodes()[j];
        Eigen::Vector3d const& Xj     = neighbour.Xi();
        Eigen::Vector4d const PXj     = polynomial(Xj);
        Eigen::RowVector4d const PXjT = PXj.transpose();
        scalar_type const Wij         = W(Xj);
        Eigen::Matrix4d const PXjPXjT = PXj * PXjT;
        // precompute moment matrix
        M_ += PXjPXjT * Wij;

        // precompute gradients of radial distance functions
        Eigen::Vector3d const gradWij = W.grad(Xj);
        gradM_[0] += PXjPXjT * gradWij(0);
        gradM_[1] += PXjPXjT * gradWij(1);
        gradM_[2] += PXjPXjT * gradWij(2);

        // precompute Aj matrices and their derivatives for meshless shape function precomputations
        // later on
        Eigen::Vector4d const Aj = PXjT * Wij;
        std::array<Eigen::Vector4d, 3u> grad_Aj{};
        grad_Aj[0] = PXjT * gradWij(0);
        grad_Aj[1] = PXjT * gradWij(1);
        grad_Aj[2] = PXjT * gradWij(2);

        Wijs_.push_back(Wij);
        gradWijs_.push_back(gradWij);
        Ajs.push_back(Aj);
        gradAjs.push_back(grad_Aj);

        // compute meshless particle density
        Vi_ += Wij;
    }

    bool is_moment_matrix_invertible{false};
    double constexpr eps = 1e-18;
    M_.computeInverseWithCheck(Minv_, is_moment_matrix_invertible, eps);
    assert(is_moment_matrix_invertible);

    Eigen::Vector3d const& Xi = W.xi();
    Eigen::Vector4d const PXi = polynomial(Xi);
    Eigen::Vector4d mesh_PXj_phi_Xis{0., 0., 0., 0.};
    std::array<Eigen::Vector4d, 3u> mesh_PXj_grad_phi_Xis{};
    mesh_PXj_grad_phi_Xis.fill(Eigen::Vector4d{0., 0., 0., 0.});
    if (body_.is_boundary_mesh_tetrahedron(ti_))
    {
        topology::tetrahedron_t const& t = body_.topology().tetrahedron(ti_);
        for (std::uint8_t i = 0u; i < 4u; ++i)
        {
            index_type const vi = t.vertex_indices()[i];

            // boundary vertices have no shape function
            if (body_.is_boundary_mesh_vertex(vi))
                continue;

            Eigen::Vector3d const& Xj = body_.x0()[vi];
            Eigen::Vector4d const PXj = polynomial(Xj);
            Eigen::Vector4d const phi = body_.phi_i(ti_, i);
            scalar_type const phi_i   = phi.dot(PXi);
            mesh_PXj_phi_Xis += PXj * phi_i;

            mesh_phi_js_[i]                  = phi_i;
            Eigen::Vector3d const grad_phi_j = body_.grad_phi_i(ti, i);
            mesh_grad_phi_js_[i]             = grad_phi_j;

            mesh_PXj_grad_phi_Xis[0] += PXj * grad_phi_j(0);
            mesh_PXj_grad_phi_Xis[1] += PXj * grad_phi_j(1);
            mesh_PXj_grad_phi_Xis[2] += PXj * grad_phi_j(2);

            // this Wij is not stored, since we're computing distance between
            // this meshless particle and a mesh particle. The mesh particle
            // is not part of the domain of the meshless shape function. However,
            // we still use the mesh particle in the density computation.
            scalar_type const Wij = W(Xj);
            Vi_ += Wij;
        }
    }
    Vi_ = static_cast<scalar_type>(1.) / Vi_;

    // Compute alpha_i, the MLS shape function coefficients
    alpha_i_ = Minv_ * (PXi - mesh_PXj_phi_Xis);

    // Compute gradient of alpha_i_
    std::array<Eigen::Vector4d, 3u> const gradPX = {
        Eigen::Vector4d{0., 1., 0., 0.},
        Eigen::Vector4d{0., 0., 1., 0.},
        Eigen::Vector4d{0., 0., 0., 1.}};
    std::array<Eigen::Vector4d, 3u> grad_Minv_PX{};
    std::array<Eigen::Vector4d, 3u> grad_Minv_sum_mesh_PXj_phi_Xis{};
    std::array<Eigen::Vector4d, 3u> grad_alpha_i{};

    for (std::uint8_t l = 0u; l < 3u; ++l)
    {
        gradMinv_[l]    = -Minv_ * gradM_[l] * Minv_;
        grad_Minv_PX[l] = (gradMinv_[l] * PXi) + (Minv_ * gradPX[l]);
        grad_Minv_sum_mesh_PXj_phi_Xis[l] =
            (gradMinv_[l] * mesh_PXj_phi_Xis) + (Minv_ * mesh_PXj_grad_phi_Xis[l]);
        grad_alpha_i[l] = grad_Minv_PX[l] - grad_Minv_sum_mesh_PXj_phi_Xis[l];
    }

    // Precompute meshless shape functions and their gradients evaluated at Xi
    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        Eigen::Vector4d const& Aj = Ajs[a];
        scalar_type const phi_j   = Aj.dot(alpha_i_);
        phi_js_.push_back(phi_j);

        Eigen::Vector3d grad_phi_j{};
        for (std::uint8_t l = 0u; l < 3u; ++l)
        {
            grad_phi_j(l) = gradAjs[a][l].dot(alpha_i_) + Aj.dot(grad_alpha_i[l]);
        }

        grad_phi_js_.push_back(grad_phi_j);
    }
}

bool hybrid_mesh_meshless_mls_node_t::is_mixed_particle() const
{
    if (ti_ == std::numeric_limits<index_type>::max())
        return false;

    bool const is_in_boundary_tetrahedron = body_.is_boundary_mesh_tetrahedron(ti_);
    return is_in_boundary_tetrahedron;
}

index_type hybrid_mesh_meshless_mls_node_t::Ni() const
{
    return ni_;
}

functions::poly6_kernel_t const& hybrid_mesh_meshless_mls_node_t::kernel() const
{
    return kernel_;
}

functions::poly6_kernel_t& hybrid_mesh_meshless_mls_node_t::kernel()
{
    return kernel_;
}

hybrid_mesh_meshless_mls_body_t const& hybrid_mesh_meshless_mls_node_t::body() const
{
    return body_;
}

std::vector<index_type> const& hybrid_mesh_meshless_mls_node_t::neighbours() const
{
    return neighbours_;
}

std::vector<scalar_type> const& hybrid_mesh_meshless_mls_node_t::Wijs() const
{
    return Wijs_;
}

std::vector<Eigen::Vector3d> const& hybrid_mesh_meshless_mls_node_t::gradWijs() const
{
    return gradWijs_;
}

scalar_type hybrid_mesh_meshless_mls_node_t::Vi() const
{
    return Vi_;
}

Eigen::Vector4d const& hybrid_mesh_meshless_mls_node_t::alphaI() const
{
    return alpha_i_;
}

std::vector<scalar_type> const& hybrid_mesh_meshless_mls_node_t::phi_js() const
{
    return phi_js_;
}

std::vector<Eigen::Vector3d> const& hybrid_mesh_meshless_mls_node_t::grad_phi_js() const
{
    return grad_phi_js_;
}

std::array<std::optional<scalar_type>, 4u> const&
hybrid_mesh_meshless_mls_node_t::mesh_phi_js() const
{
    return mesh_phi_js_;
}

std::array<std::optional<Eigen::Vector3d>, 4u> const&
hybrid_mesh_meshless_mls_node_t::mesh_grad_phi_js() const
{
    return mesh_grad_phi_js_;
}

Eigen::Vector3d const& hybrid_mesh_meshless_mls_node_t::Xi() const
{
    return kernel_.xi();
}

Eigen::Vector3d& hybrid_mesh_meshless_mls_node_t::Xi()
{
    return kernel_.xi();
}

Eigen::Matrix3d const& hybrid_mesh_meshless_mls_node_t::Fi() const
{
    return Fi_;
}

Eigen::Matrix3d& hybrid_mesh_meshless_mls_node_t::Fi()
{
    return Fi_;
}

Eigen::Vector3d const& hybrid_mesh_meshless_mls_node_t::xi() const
{
    return xi_;
}

Eigen::Vector3d& hybrid_mesh_meshless_mls_node_t::xi()
{
    return xi_;
}

Eigen::Matrix3d const& hybrid_mesh_meshless_mls_node_t::Ei() const
{
    return Ei_;
}

Eigen::Matrix3d& hybrid_mesh_meshless_mls_node_t::Ei()
{
    return Ei_;
}

index_type hybrid_mesh_meshless_mls_node_t::ti() const
{
    return ti_;
}

std::array<index_type, 4u> const& hybrid_mesh_meshless_mls_node_t::vis() const
{
    topology::tetrahedron_t const& t = body_.topology().tetrahedron(ti_);
    return t.vertex_indices();
}

Eigen::Vector4d hybrid_mesh_meshless_mls_node_t::polynomial(Eigen::Vector3d const& X) const
{
    return Eigen::Vector4d(1., X.x(), X.y(), X.z());
}

} // namespace mechanics
} // namespace physics
} // namespace sbs