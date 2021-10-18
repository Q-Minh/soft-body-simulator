#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_body.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_node.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_surface.h>

namespace sbs {
namespace physics {
namespace mechanics {

hybrid_mesh_meshless_sph_surface_t::hybrid_mesh_meshless_sph_surface_t(
    hybrid_mesh_meshless_sph_body_t* mechanical_model)
    : tetrahedral_mesh_boundary_t(&mechanical_model->topology()),
      world_space_vertices_(vertex_count()),
      Xkjs_(),
      Wkjs_(),
      sks_(),
      neighbours_(),
      mechanical_model_(mechanical_model)
{
}

void hybrid_mesh_meshless_sph_surface_t::initialize_interpolation_scheme()
{
    Xkjs_.clear();
    Wkjs_.clear();
    sks_.clear();
    neighbours_.clear();
    Xkjs_.resize(vertex_count());
    Wkjs_.resize(vertex_count());
    sks_.resize(vertex_count());
    neighbours_.resize(vertex_count());

    detail::hybrid_mesh_meshless_sph::meshless_node_range_searcher_t const&
        meshless_range_searcher = mechanical_model_->meshless_node_range_searcher();
    std::vector<hybrid_mesh_meshless_sph_node_t> const& nodes = mechanical_model_->meshless_nodes();
    scalar_type const h                                       = mechanical_model_->h();

    for (std::size_t i = 0u; i < vertex_count(); ++i)
    {
        Eigen::Vector3d const& Xk = material_space_vertex(i).position;
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

            Xkjs_[i].push_back(Xkj);
            Wkjs_[i].push_back(Wkj);
            neighbours_[i].push_back(j);

            sk += Vj * Wkj;
        }
        sk      = 1. / sk;
        sks_[i] = sk;

        world_space_vertices_[i].position = material_space_vertex(i).position;
        world_space_vertices_[i].color    = material_space_vertex(i).color;
    }
}

hybrid_mesh_meshless_sph_surface_t::vertex_type
hybrid_mesh_meshless_sph_surface_t::vertex(std::size_t vi) const
{
    return world_space_vertices_[vi];
}

void hybrid_mesh_meshless_sph_surface_t::prepare_vertices_for_rendering()
{
    prepare_vertices_for_surface_rendering();
}

hybrid_mesh_meshless_sph_surface_t::vertex_type&
hybrid_mesh_meshless_sph_surface_t::world_space_vertex(std::size_t vi)
{
    return world_space_vertices_[vi];
}

hybrid_mesh_meshless_sph_surface_t::vertex_type&
hybrid_mesh_meshless_sph_surface_t::material_space_vertex(std::size_t vi)
{
    return mutable_vertex(vi);
}

void hybrid_mesh_meshless_sph_surface_t::compute_positions()
{
    /**
     * IMPORTANT NOTE:
     * Because the surface vertices correspond exactly to boundary vertices of
     * the initial tetrahedral mesh, then the mesh nodes' shape functions all evaluate
     * to 0 at these surface vertices (remember, no boundary vertex has a shape function
     * in the hybrid model). Because of this, the interpolation to compute positions is
     * not affected by the mesh nodes' shape functions that overlap with the mixed domain
     * of particles and mesh nodes. This is why we only compute positions using neighboring
     * meshless nodes. However, if we choose this embedded surface to be any other way,
     * we would have to also add the contributions of the mesh nodes' shape functions
     * to this interpolation to find surface mesh vertex positions. This should be done
     * in the future, since we want to be able to embed high resolution surfaces into
     * this hybrid mesh/meshless model instead of always using the initial tetrahedral
     * mesh's boundary.
     */
    auto const& meshless_nodes = mechanical_model_->meshless_nodes();
    for (std::size_t i = 0u; i < vertex_count(); ++i)
    {
        Eigen::Vector3d& xk = world_space_vertices_[i].position;
        xk.setZero();

        auto const num_neighbours = Xkjs_[i].size();
        for (std::size_t b = 0u; b < num_neighbours; ++b)
        {
            index_type const j                        = neighbours_[i][b];
            hybrid_mesh_meshless_sph_node_t const& nj = meshless_nodes[j];
            Eigen::Vector3d const& Xkj                = Xkjs_[i][b];
            scalar_type const Wkj                     = Wkjs_[i][b];
            scalar_type const Vj                      = nj.Vi();
            Eigen::Matrix3d const& Fj                 = nj.Fi();
            Eigen::Vector3d const& xj                 = nj.xi();
            xk += Vj * (Fj * Xkj + xj) * Wkj;
        }

        scalar_type const sk = sks_[i];
        xk                   = sk * xk;
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

        Eigen::Vector3d const& p1 = world_space_vertices_[v1].position;
        Eigen::Vector3d const& p2 = world_space_vertices_[v2].position;
        Eigen::Vector3d const& p3 = world_space_vertices_[v3].position;

        Eigen::Vector3d const n = (p2 - p1).cross(p3 - p1);

        world_space_vertices_[v1].normal += n;
        world_space_vertices_[v2].normal += n;
        world_space_vertices_[v3].normal += n;
    }

    for (std::size_t i = 0u; i < world_space_vertices_.size(); ++i)
    {
        world_space_vertices_[i].normal.normalize();
    }
}

void hybrid_mesh_meshless_sph_surface_t::prepare_vertices_for_surface_rendering()
{
    std::size_t constexpr num_attributes_per_vertex = 9u;
    std::size_t const vertex_count                  = world_space_vertices_.size();
    std::vector<float> vertex_buffer{};
    vertex_buffer.reserve(vertex_count * num_attributes_per_vertex);

    for (vertex_type const& vertex : world_space_vertices_)
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
