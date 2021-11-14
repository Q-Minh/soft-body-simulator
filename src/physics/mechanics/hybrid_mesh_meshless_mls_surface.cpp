#include <sbs/physics/mechanics/hybrid_mesh_meshless_mls_body.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_mls_node.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_mls_surface.h>

namespace sbs {
namespace physics {
namespace mechanics {

hybrid_mesh_meshless_sph_surface_t::hybrid_mesh_meshless_sph_surface_t(
    hybrid_mesh_meshless_mls_body_t* mechanical_model,
    std::vector<Eigen::Vector3d> const& vertices,
    std::vector<topology::triangle_t> const& triangles)
    : render_vertices_(), triangles_(), vertices_(), mechanical_model_(mechanical_model)
{
    render_vertices_.reserve(vertices.size());
    vertices_.reserve(vertices.size());
    for (auto const& p : vertices)
    {
        vertex_type v{};
        v.position = p;
        render_vertices_.push_back(v);
        hybrid_mesh_meshless_mls_surface_vertex_t mv(p);
        vertices_.push_back(mv);
    }
    triangles_.reserve(triangles.size());
    for (auto const& triangle : triangles)
    {
        triangle_type f{};
        std::copy(
            triangle.vertex_indices().begin(),
            triangle.vertex_indices().end(),
            f.vertices.begin());

        triangles_.push_back(f);
    }
}

void hybrid_mesh_meshless_sph_surface_t::initialize_interpolation_scheme(scalar_type const h)
{
    detail::hybrid_mesh_meshless_sph::meshless_node_range_searcher_t const&
        meshless_range_searcher = mechanical_model_->meshless_node_range_searcher();
    detail::hybrid_mesh_meshless_sph::mesh_tetrahedron_range_searcher_t const&
        mesh_tet_range_searcher = mechanical_model_->mesh_tetrahedron_range_searcher();

    std::vector<hybrid_mesh_meshless_mls_node_t> const& nodes = mechanical_model_->meshless_nodes();

    for (std::size_t i = 0u; i < vertices_.size(); ++i)
    {
        Eigen::Vector3d const& Xk = vertices_[i].Xi();
        std::vector<index_type> const& neighbour_indices =
            meshless_range_searcher.neighbours_of(Xk, h);

        index_type const ti = mesh_tet_range_searcher.in_tetrahedron(Xk);
        vertices_[i].initialize(neighbour_indices, *mechanical_model_, ti);
        render_vertices_[i].position = Xk;
    }
}

std::size_t hybrid_mesh_meshless_sph_surface_t::vertex_count() const
{
    return render_vertices_.size();
}

std::size_t hybrid_mesh_meshless_sph_surface_t::triangle_count() const
{
    return triangles_.size();
}

hybrid_mesh_meshless_sph_surface_t::vertex_type
hybrid_mesh_meshless_sph_surface_t::vertex(std::size_t vi) const
{
    return render_vertices_[vi];
}

common::shared_vertex_surface_mesh_i::triangle_type
hybrid_mesh_meshless_sph_surface_t::triangle(std::size_t f) const
{
    return triangles_[f];
}

void hybrid_mesh_meshless_sph_surface_t::prepare_vertices_for_rendering()
{
    prepare_vertices_for_surface_rendering();
}

void hybrid_mesh_meshless_sph_surface_t::prepare_indices_for_rendering()
{
    std::size_t const index_count = triangles_.size() * 3u;
    std::vector<std::uint32_t> index_buffer{};
    index_buffer.reserve(index_count);
    for (triangle_type const& f : triangles_)
    {
        index_buffer.push_back(f.vertices[0]);
        index_buffer.push_back(f.vertices[1]);
        index_buffer.push_back(f.vertices[2]);
    }

    transfer_indices_for_rendering(std::move(index_buffer));
}

hybrid_mesh_meshless_sph_surface_t::vertex_type&
hybrid_mesh_meshless_sph_surface_t::world_space_vertex(std::size_t vi)
{
    return render_vertices_[vi];
}

Eigen::Vector3d& hybrid_mesh_meshless_sph_surface_t::material_space_position(std::size_t vi)
{
    return vertices_[vi].Xi();
}

std::vector<hybrid_mesh_meshless_mls_surface_vertex_t> const&
hybrid_mesh_meshless_sph_surface_t::embedded_vertices() const
{
    return vertices_;
}

hybrid_mesh_meshless_mls_body_t const* hybrid_mesh_meshless_sph_surface_t::mechanical_model() const
{
    return mechanical_model_;
}

hybrid_mesh_meshless_mls_body_t* hybrid_mesh_meshless_sph_surface_t::mechanical_model()
{
    return mechanical_model_;
}

void hybrid_mesh_meshless_sph_surface_t::compute_positions()
{
    auto const& nodes = mechanical_model_->meshless_nodes();
    for (std::size_t i = 0u; i < vertex_count(); ++i)
    {
        Eigen::Vector3d xk{0., 0., 0.};

        auto& meshless_surface_vertex = vertices_[i];
        auto const& phi_js            = meshless_surface_vertex.phi_js();

        auto const& neighbours = meshless_surface_vertex.neighbours();
        for (std::size_t b = 0u; b < neighbours.size(); ++b)
        {
            index_type const j                        = neighbours[b];
            hybrid_mesh_meshless_mls_node_t const& nj = nodes[j];
            Eigen::Vector3d const& xj                 = nj.xi();
            scalar_type const phi_j                   = phi_js[b];
            xk += xj * phi_j;
        }
        if (meshless_surface_vertex.is_in_tetrahedron())
        {
            auto const& mesh_phi_js = meshless_surface_vertex.mesh_phi_js();
            topology::tetrahedron_t const& t =
                mechanical_model_->topology().tetrahedron(meshless_surface_vertex.ti());
            auto const& vis = t.vertex_indices();
            for (std::uint8_t v = 0u; v < 4u; ++v)
            {
                if (!mesh_phi_js[v].has_value())
                    continue;

                index_type const vi       = vis[v];
                scalar_type const phi_j   = mesh_phi_js[v].value();
                Eigen::Vector3d const& xj = mechanical_model_->x()[vi];
                xk += xj * phi_j;
            }
        }

        meshless_surface_vertex.xi() = xk;
        render_vertices_[i].position = xk;
    }
}

void hybrid_mesh_meshless_sph_surface_t::compute_normals()
{
    for (std::size_t i = 0u; i < triangle_count(); ++i)
    {
        triangle_type const& f = triangle(i);
        auto const v1          = f.vertices[0u];
        auto const v2          = f.vertices[1u];
        auto const v3          = f.vertices[2u];

        Eigen::Vector3d const& p1 = render_vertices_[v1].position;
        Eigen::Vector3d const& p2 = render_vertices_[v2].position;
        Eigen::Vector3d const& p3 = render_vertices_[v3].position;

        Eigen::Vector3d const n = (p2 - p1).cross(p3 - p1);

        render_vertices_[v1].normal += n;
        render_vertices_[v2].normal += n;
        render_vertices_[v3].normal += n;
    }

    for (std::size_t i = 0u; i < render_vertices_.size(); ++i)
    {
        render_vertices_[i].normal.normalize();
    }
}

void hybrid_mesh_meshless_sph_surface_t::prepare_vertices_for_surface_rendering()
{
    std::size_t constexpr num_attributes_per_vertex = 9u;
    std::size_t const vertex_count                  = render_vertices_.size();
    std::vector<float> vertex_buffer{};
    vertex_buffer.reserve(vertex_count * num_attributes_per_vertex);

    for (vertex_type const& vertex : render_vertices_)
    {
        vertex_buffer.push_back(static_cast<float>(vertex.position.x()));
        vertex_buffer.push_back(static_cast<float>(vertex.position.y()));
        vertex_buffer.push_back(static_cast<float>(vertex.position.z()));

        // use triangle normal since face-based
        vertex_buffer.push_back(static_cast<float>(vertex.normal.x()));
        vertex_buffer.push_back(static_cast<float>(vertex.normal.y()));
        vertex_buffer.push_back(static_cast<float>(vertex.normal.z()));

        vertex_buffer.push_back(vertex.color.x());
        vertex_buffer.push_back(vertex.color.y());
        vertex_buffer.push_back(vertex.color.z());
    }

    transfer_vertices_for_rendering(std::move(vertex_buffer));
}

void hybrid_mesh_meshless_mls_surface_vertex_t::initialize(
    std::vector<index_type> const& meshless_neighbours,
    hybrid_mesh_meshless_mls_body_t const& mechanical_model,
    index_type ti)
{
    ti_ = ti;
    bool const is_boundary_tet =
        is_in_tetrahedron() && mechanical_model.is_boundary_mesh_tetrahedron(ti_);

    neighbours_ = meshless_neighbours;

    auto const& Xi = Xi_;
    Eigen::Vector4d const PXi{1., Xi.x(), Xi.y(), Xi.z()};
    phi_js_.clear();
    mesh_phi_js_.fill({});

    // First, compute moment matrix and our Ajs
    Eigen::Matrix4d M{};
    M.setZero();
    std::vector<Eigen::RowVector4d> Ajs{};
    Ajs.reserve(meshless_neighbours.size());
    for (std::size_t a = 0u; a < meshless_neighbours.size(); ++a)
    {
        index_type const j                               = meshless_neighbours[a];
        hybrid_mesh_meshless_mls_node_t const& neighbour = mechanical_model.meshless_nodes()[j];
        Eigen::Vector3d const& Xj                        = neighbour.Xi();
        auto const& W                                    = neighbour.kernel();
        Eigen::Vector4d const PXj                        = neighbour.polynomial(Xj);
        Eigen::RowVector4d const PXjT                    = PXj.transpose();
        scalar_type const Wij                            = W(Xj);
        Eigen::Matrix4d const Minc                       = PXj * PXjT * Wij;
        M += Minc;

        Eigen::RowVector4d const Aj = PXjT * Wij;
        Ajs.push_back(Aj);
    }

    bool is_moment_matrix_invertible{false};
    double constexpr eps = 1e-18;
    Eigen::Matrix4d Minv{};
    M.computeInverseWithCheck(Minv, is_moment_matrix_invertible, eps);
    assert(is_moment_matrix_invertible);

    // initialize rhs for solving for alpha
    Eigen::Vector4d b{1., Xi.x(), Xi.y(), Xi.z()};

    // update rhs if there are mesh shape functions
    if (is_boundary_tet)
    {
        topology::tetrahedron_t const& t = mechanical_model.topology().tetrahedron(ti_);
        for (std::uint8_t i = 0u; i < 4u; ++i)
        {
            index_type const vi = t.vertex_indices()[i];
            // boundary vertices have no shape function
            if (mechanical_model.is_boundary_mesh_vertex(vi))
                continue;

            Eigen::Vector3d const& Xj = mechanical_model.x0()[vi];
            Eigen::Vector4d const PXj{1., Xj.x(), Xj.y(), Xj.z()};
            Eigen::Vector4d const& phi_i = mechanical_model.phi_i(ti_, i);
            scalar_type const mesh_phi_j = phi_i.dot(PXi);
            b -= PXj * mesh_phi_j;

            // store these mesh shape functions for later use
            mesh_phi_js_[i] = mesh_phi_j;
        }
    }

    Eigen::Vector4d const alpha = Minv * b;

    // Precompute meshless shape functions
    for (std::size_t a = 0u; a < meshless_neighbours.size(); ++a)
    {
        Eigen::RowVector4d const& Aj = Ajs[a];
        scalar_type const phi_j      = Aj.dot(alpha);
        phi_js_.push_back(phi_j);
    }
}

Eigen::Vector3d const& hybrid_mesh_meshless_mls_surface_vertex_t::xi() const
{
    return xi_;
}

Eigen::Vector3d& hybrid_mesh_meshless_mls_surface_vertex_t::xi()
{
    return xi_;
}

Eigen::Vector3d const& hybrid_mesh_meshless_mls_surface_vertex_t::Xi() const
{
    return Xi_;
}

Eigen::Vector3d& hybrid_mesh_meshless_mls_surface_vertex_t::Xi()
{
    return Xi_;
}

std::vector<scalar_type> const& hybrid_mesh_meshless_mls_surface_vertex_t::phi_js() const
{
    return phi_js_;
}

std::array<std::optional<scalar_type>, 4u> const&
hybrid_mesh_meshless_mls_surface_vertex_t::mesh_phi_js() const
{
    return mesh_phi_js_;
}

std::vector<index_type> const& hybrid_mesh_meshless_mls_surface_vertex_t::neighbours() const
{
    return neighbours_;
}

} // namespace mechanics
} // namespace physics
} // namespace sbs
