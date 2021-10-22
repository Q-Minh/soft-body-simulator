#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_body.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_node.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_surface.h>

namespace sbs {
namespace physics {
namespace mechanics {

hybrid_mesh_meshless_sph_surface_t::hybrid_mesh_meshless_sph_surface_t(
    hybrid_mesh_meshless_sph_body_t* mechanical_model,
    std::vector<Eigen::Vector3d> const& vertices,
    std::vector<triangle_t> const& triangles)
    : render_vertices_(), triangles_(), vertices_(), mechanical_model_(mechanical_model)
{
    render_vertices_.reserve(vertices.size());
    vertices_.reserve(vertices.size());
    for (auto const& p : vertices)
    {
        vertex_type v{};
        v.position = p;
        render_vertices_.push_back(v);
        hybrid_mesh_meshless_sph_surface_vertex_t mv(p);
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

    std::vector<hybrid_mesh_meshless_sph_node_t> const& nodes = mechanical_model_->meshless_nodes();

    for (std::size_t i = 0u; i < vertices_.size(); ++i)
    {
        std::vector<Eigen::Vector3d> Xkjs{};
        std::vector<scalar_type> Wkjs{};
        std::vector<scalar_type> Vjs{};
        std::vector<index_type> neighbours{};

        Eigen::Vector3d const& Xk = vertices_[i].x0();
        std::vector<index_type> const& neighbour_indices =
            meshless_range_searcher.neighbours_of(Xk, h);

        scalar_type sk{0.};
        for (std::size_t a = 0u; a < neighbour_indices.size(); ++a)
        {
            index_type const j            = neighbour_indices[a];
            meshless_sph_node_t const& nj = nodes[j];
            scalar_type const Vj          = nj.Vi();
            Eigen::Vector3d const Xj      = nj.Xi();
            Eigen::Vector3d const& Xkj    = Xk - Xj;
            scalar_type const Wkj         = nj.kernel()(Xk);

            Xkjs.push_back(Xkj);
            Wkjs.push_back(Wkj);
            Vjs.push_back(Vj);
            neighbours.push_back(j);

            sk += Vj * Wkj;
        }
        sk = 1. / sk;

        index_type const ti = mesh_tet_range_searcher.in_tetrahedron(Xk);
        vertices_[i] =
            hybrid_mesh_meshless_sph_surface_vertex_t(Xk, Xk, Xkjs, Wkjs, Vjs, sk, neighbours, ti);
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
    return vertices_[vi].x0();
}

std::vector<hybrid_mesh_meshless_sph_surface_vertex_t> const&
hybrid_mesh_meshless_sph_surface_t::embedded_vertices() const
{
    return vertices_;
}

hybrid_mesh_meshless_sph_body_t const* hybrid_mesh_meshless_sph_surface_t::mechanical_model() const
{
    return mechanical_model_;
}

hybrid_mesh_meshless_sph_body_t* hybrid_mesh_meshless_sph_surface_t::mechanical_model()
{
    return mechanical_model_;
}

void hybrid_mesh_meshless_sph_surface_t::compute_positions()
{
    auto const& nodes = mechanical_model_->meshless_nodes();
    for (std::size_t i = 0u; i < vertex_count(); ++i)
    {
        Eigen::Vector3d& xk = render_vertices_[i].position;
        xk.setZero();

        auto const& neighbours = vertices_[i].neighbours();
        auto const& Xkjs       = vertices_[i].Xkjs();
        auto const& Wkjs       = vertices_[i].Wkjs();
        auto const& Vjs        = vertices_[i].Vjs();
        auto const& sk         = vertices_[i].sk();
        for (std::size_t b = 0u; b < neighbours.size(); ++b)
        {
            index_type const j            = neighbours[b];
            meshless_sph_node_t const& nj = nodes[j];
            Eigen::Vector3d const& Xkj    = Xkjs[b];
            scalar_type const Wkj         = Wkjs[b];
            scalar_type const Vj          = Vjs[b];
            Eigen::Matrix3d const& Fj     = nj.Fi();
            Eigen::Vector3d const& xj     = nj.xi();
            xk += Vj * (Fj * Xkj + xj) * Wkj;
        }
        xk               = sk * xk;
        vertices_[i].x() = xk;

        // auto const ti                         = meshless_vertex.ti();
        // bool const is_vertex_in_boundary_mesh = ti != std::numeric_limits<index_type>::max() &&
        //                                         mechanical_model_->is_boundary_mesh_tetrahedron(ti);
        // if (is_vertex_in_boundary_mesh)
        //{
        //     // Compute contribution of interpolation over mesh nodes
        //     tetrahedron_t const& t = mechanical_model_->topology().tetrahedron(ti);
        //     for (std::uint8_t v = 0u; v < 4u; ++v)
        //     {
        //         index_type const vi = t.vertex_indices()[v];
        //         if (mechanical_model_->is_boundary_mesh_vertex(vi))
        //             continue;

        //        auto const& xi   = mechanical_model_->x()[vi];
        //        auto const phi   = mechanical_model_->phi_i(ti, v);
        //        auto const& Xi   = meshless_vertex.x0();
        //        auto const phi_i = phi.dot(Eigen::Vector4d{1., Xi.x(), Xi.y(), Xi.z()});
        //        xk += xi * phi_i;
        //    }
        //}
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

} // namespace mechanics
} // namespace physics
} // namespace sbs
