#include <Discregrid/cubic_lagrange_discrete_grid.hpp>
#include <Discregrid/geometry/mesh_distance.hpp>
#include <Discregrid/mesh/triangle_mesh.hpp>
#include <algorithm>
#include <sbs/common/geometry.h>
#include <sbs/physics/mechanics/meshless_sph_body.h>
#include <sbs/physics/mechanics/meshless_sph_node.h>
#include <sbs/physics/particle.h>
#include <sbs/physics/simulation.h>

namespace sbs {
namespace physics {
namespace mechanics {

meshless_sph_body_t::meshless_sph_body_t(
    simulation_t& simulation,
    index_type id,
    common::geometry_t const& geometry,
    scalar_type const h,
    std::array<unsigned int, 3u> const& resolution)
    : body_t(simulation, id),
      meshless_nodes_(),
      material_space_range_query_(),
      volumetric_topology_(),
      visual_model_(),
      collision_model_(),
      h_()
{
    assert(geometry.geometry_type == common::geometry_t::geometry_type_t::tetrahedron);
    assert(geometry.has_indices());
    assert(geometry.has_positions());

    /**
     * Compute topology of the given tetrahedral mesh geometry
     */
    volumetric_topology_.reserve_vertices(geometry.positions.size() / 3u);

    for (std::size_t i = 0u; i < geometry.indices.size(); i += 4u)
    {
        sbs::physics::tetrahedron_t const tetrahedron{
            static_cast<index_type>(geometry.indices[i]),
            static_cast<index_type>(geometry.indices[i + 1u]),
            static_cast<index_type>(geometry.indices[i + 2u]),
            static_cast<index_type>(geometry.indices[i + 3u])};

        volumetric_topology_.add_tetrahedron(tetrahedron);
    }
    /**
     * Compute boundary surface mesh of the given tetrahedral mesh geometry
     */
    std::vector<Eigen::Vector3d> vertices{};
    vertices.reserve(volumetric_topology_.vertex_count());
    for (std::size_t i = 0u; i < volumetric_topology_.vertex_count(); ++i)
    {
        auto const idx      = i * 3u;
        scalar_type const x = static_cast<scalar_type>(geometry.positions[idx]);
        scalar_type const y = static_cast<scalar_type>(geometry.positions[idx + 1u]);
        scalar_type const z = static_cast<scalar_type>(geometry.positions[idx + 2u]);
        Eigen::Vector3d const pos{x, y, z};
        vertices.push_back(pos);
    }
    std::vector<std::array<unsigned int, 3u>> faces{};
    std::vector<triangle_t> const boundary_triangles =
        volumetric_topology_.oriented_boundary_triangles();
    faces.reserve(boundary_triangles.size());
    for (triangle_t const& f : boundary_triangles)
        faces.push_back({f.v1(), f.v2(), f.v3()});

    /**
     * Compute the SDF of that boundary surface to obtain an inside/outside predicate
     * for our particle sampling
     */
    Discregrid::TriangleMesh mesh(vertices, faces);
    Eigen::AlignedBox3d domain{};
    for (auto const& x : mesh.vertices())
    {
        for (auto const& x : mesh.vertices())
        {
            domain.extend(x);
        }
    }

    // add a bit of tolerance to the domain to fully enclose the mesh
    domain.min() -= Eigen::Vector3d{1e-3, 1e-3, 1e-3};
    domain.max() += Eigen::Vector3d{1e-3, 1e-3, 1e-3};

    Discregrid::MeshDistance md(mesh);

    auto const sdf = [&md](Eigen::Vector3d const& xi) {
        return md.signedDistanceCached(xi);
    };

    Discregrid::CubicLagrangeDiscreteGrid grid(domain, resolution);
    auto const sdf_id = grid.addFunction(sdf);

    // Compute grid cell dimensions for the grid enclosing the boundary surface mesh
    scalar_type const dx =
        (domain.max().x() - domain.min().x()) / static_cast<scalar_type>(resolution[0]);
    scalar_type const dy =
        (domain.max().y() - domain.min().y()) / static_cast<scalar_type>(resolution[1]);
    scalar_type const dz =
        (domain.max().z() - domain.min().z()) / static_cast<scalar_type>(resolution[2]);

    // The support radius of the meshless meshless_nodes is a multiple h times the largest grid cell
    // dimension
    h_ = h * std::max({dx, dy, dz});
    h_ += std::numeric_limits<scalar_type>::epsilon();

    // Create a meshless node/particle for each grid cell center position
    // that lies inside the boundary surface using its SDF
    for (unsigned int i = 0u; i < resolution[0]; ++i)
    {
        for (unsigned int j = 0u; j < resolution[1]; ++j)
        {
            for (unsigned int k = 0u; k < resolution[2]; ++k)
            {
                // Place grid node at the center of grid cells
                scalar_type const x =
                    domain.min().x() + static_cast<scalar_type>(i) * dx + 0.5 * dx;
                scalar_type const y =
                    domain.min().y() + static_cast<scalar_type>(j) * dy + 0.5 * dy;
                scalar_type const z =
                    domain.min().z() + static_cast<scalar_type>(k) * dz + 0.5 * dz;

                Eigen::Vector3d const& grid_node_position{x, y, z};
                scalar_type const signed_distance = grid.interpolate(sdf_id, grid_node_position);
                bool const is_grid_node_inside_surface = signed_distance < 0.;
                if (!is_grid_node_inside_surface)
                    continue;

                functions::poly6_kernel_t kernel(grid_node_position, h_);
                auto const ni = static_cast<index_type>(meshless_nodes_.size());
                meshless_sph_node_t node(ni, kernel);
                node.xi() = kernel.xi();
                meshless_nodes_.push_back(node);
            }
        }
    }

    // Compress the mesh (vertices, faces) into a more compact mesh (surface_vertices,
    // surface_triangles) where the surface_vertices vector only contains elements in vertices that
    // are actually referenced by faces
    std::unordered_map<index_type, index_type> tet_vertex_to_surface_vertex{};
    std::vector<Eigen::Vector3d> surface_vertices{};
    std::vector<triangle_t> surface_triangles{};
    for (auto const& f : boundary_triangles)
    {
        if (tet_vertex_to_surface_vertex.find(f.v1()) == tet_vertex_to_surface_vertex.end())
        {
            auto const index                     = tet_vertex_to_surface_vertex.size();
            tet_vertex_to_surface_vertex[f.v1()] = static_cast<index_type>(index);
            surface_vertices.push_back(vertices[f.v1()]);
        }
        if (tet_vertex_to_surface_vertex.find(f.v2()) == tet_vertex_to_surface_vertex.end())
        {
            auto const index                     = tet_vertex_to_surface_vertex.size();
            tet_vertex_to_surface_vertex[f.v2()] = static_cast<index_type>(index);
            surface_vertices.push_back(vertices[f.v2()]);
        }
        if (tet_vertex_to_surface_vertex.find(f.v3()) == tet_vertex_to_surface_vertex.end())
        {
            auto const index                     = tet_vertex_to_surface_vertex.size();
            tet_vertex_to_surface_vertex[f.v3()] = static_cast<index_type>(index);
            surface_vertices.push_back(vertices[f.v3()]);
        }

        auto const v1 = tet_vertex_to_surface_vertex[f.v1()];
        auto const v2 = tet_vertex_to_surface_vertex[f.v2()];
        auto const v3 = tet_vertex_to_surface_vertex[f.v3()];
        triangle_t const triangle{v1, v2, v3};
        surface_triangles.push_back(triangle);
    }

    // Compute the initial visual mesh as the boundary surface mesh of the given geometry
    // and assign colors to its vertices
    visual_model_ = meshless_sph_surface_t(this, surface_vertices, surface_triangles);
    for (std::size_t i = 0u; i < vertices.size(); ++i)
    {
        if (tet_vertex_to_surface_vertex.find(static_cast<index_type>(i)) ==
            tet_vertex_to_surface_vertex.end())
            continue;

        auto const idx = i * 3u;
        float const r  = static_cast<float>(geometry.colors[idx] / 255.f);
        float const g  = static_cast<float>(geometry.colors[idx + 1u] / 255.f);
        float const b  = static_cast<float>(geometry.colors[idx + 2u] / 255.f);

        auto const vi = tet_vertex_to_surface_vertex[static_cast<index_type>(i)];
        visual_model_.world_space_vertex(vi).color = Eigen::Vector3f{r, g, b};
    }
}

body_t::visual_model_type const& meshless_sph_body_t::visual_model() const
{
    return visual_model_;
}

body_t::collision_model_type const& meshless_sph_body_t::collision_model() const
{
    return collision_model_;
}

body_t::visual_model_type& meshless_sph_body_t::visual_model()
{
    return visual_model_;
}

body_t::collision_model_type& meshless_sph_body_t::collision_model()
{
    return collision_model_;
}

void meshless_sph_body_t::update_visual_model()
{
    visual_model_.compute_positions();
    visual_model_.compute_normals();
}

void meshless_sph_body_t::update_collision_model()
{
    collision_model_.update(simulation());
}

void meshless_sph_body_t::update_physical_model()
{
    auto const& particles = simulation().particles()[id()];
    assert(particles.size() == meshless_nodes_.size());
    for (std::size_t i = 0u; i < particles.size(); ++i)
    {
        meshless_sph_node_t& node = meshless_nodes_[i];
        particle_t const& p       = particles[i];
        node.xi()                 = p.x();
    }
}

void meshless_sph_body_t::transform(Eigen::Affine3d const& affine)
{
    for (auto& meshless_node : meshless_nodes_)
    {
        meshless_node.Xi() = affine * meshless_node.Xi().homogeneous();
        meshless_node.xi() = affine * meshless_node.xi().homogeneous();
    }

    for (std::size_t i = 0u; i < visual_model_.vertex_count(); ++i)
    {
        visual_model_.world_space_vertex(i).position =
            affine * visual_model_.world_space_vertex(i).position.homogeneous();
        visual_model_.material_space_position(i) =
            affine * visual_model_.material_space_position(i).homogeneous();
    }

    auto const v             = Eigen::Vector3d{1., 1., 1.}.normalized();
    auto const vp            = affine * Eigen::Vector4d{v.x(), v.y(), v.z(), 0.};
    auto const length_change = vp.norm();
    h_ *= length_change;
}

std::vector<meshless_sph_node_t> const& meshless_sph_body_t::nodes() const
{
    return meshless_nodes_;
}

std::vector<meshless_sph_node_t>& meshless_sph_body_t::nodes()
{
    return meshless_nodes_;
}

meshless_sph_surface_t const& meshless_sph_body_t::surface_mesh() const
{
    return visual_model_;
}

meshless_sph_surface_t& meshless_sph_body_t::surface_mesh()
{
    return visual_model_;
}

collision::point_bvh_model_t const& meshless_sph_body_t::bvh() const
{
    return collision_model_;
}

meshless_sph_body_range_searcher_t const& meshless_sph_body_t::range_searcher() const
{
    return material_space_range_query_;
}

tetrahedron_set_t const& meshless_sph_body_t::topology() const
{
    return volumetric_topology_;
}

tetrahedron_set_t& meshless_sph_body_t::topology()
{
    return volumetric_topology_;
}

scalar_type meshless_sph_body_t::h() const
{
    return h_;
}

void meshless_sph_body_t::initialize_physical_model()
{
    // Get a spatial acceleration query object to obtain neighbours of
    // each meshless node in material space
    material_space_range_query_ = meshless_sph_body_range_searcher_t(&meshless_nodes_);
    std::vector<std::vector<Eigen::Vector3d const*>> Xjs{};
    std::vector<std::vector<index_type>> node_neighbour_indices{};
    std::vector<std::vector<meshless_sph_node_t const*>> node_neighbours{};

    Xjs.resize(meshless_nodes_.size());
    node_neighbour_indices.resize(meshless_nodes_.size());
    node_neighbours.resize(meshless_nodes_.size());

    // precompute all quantities that depend only on material space
    for (std::size_t i = 0u; i < meshless_nodes_.size(); ++i)
    {
        meshless_sph_node_t const& node = meshless_nodes_[i];
        auto const ni                   = static_cast<index_type>(i);
        std::vector<index_type> const neighbours_in_domain =
            material_space_range_query_.neighbours_of(ni);

        // precompute all neighbour information required to initialize our meshless nodes
        for (std::size_t k = 0u; k < neighbours_in_domain.size(); ++k)
        {
            index_type const j                   = neighbours_in_domain[k];
            meshless_sph_node_t const& neighbour = meshless_nodes_[j];
            // neighbour positions
            Eigen::Vector3d const* Xj = &neighbour.Xi();
            Xjs[i].push_back(Xj);
            // neighbour indices
            node_neighbour_indices[i].push_back(j);
            // neighbour meshless nodes
            node_neighbours[i].push_back(&neighbour);
        }
    }
    // update meshless nodes using precomputed neighbour information
    for (std::size_t i = 0u; i < meshless_nodes_.size(); ++i)
    {
        auto const& xi = meshless_nodes_[i].xi();
        meshless_nodes_[i].initialize(xi, Xjs[i], node_neighbour_indices[i]);
    }
    // precompute the correction matrix Li for each meshless node
    for (std::size_t i = 0u; i < meshless_nodes_.size(); ++i)
    {
        meshless_nodes_[i].cache_Li_Vj(node_neighbours[i]);
    }
}

void meshless_sph_body_t::initialize_visual_model()
{
    visual_model_.initialize_interpolation_scheme(h_ * 2.);
    visual_model_.compute_normals();
}

void meshless_sph_body_t::initialize_collision_model()
{
    collision_model_      = collision::point_bvh_model_t(&visual_model_);
    collision_model_.id() = id();
}

meshless_sph_body_range_searcher_t::meshless_sph_body_range_searcher_t() : base_type(0u), nodes_()
{
}

meshless_sph_body_range_searcher_t::meshless_sph_body_range_searcher_t(
    std::vector<meshless_sph_node_t> const* nodes)
    : base_type(nodes->size()), nodes_(nodes)
{
    this->construct();
}

std::vector<index_type> meshless_sph_body_range_searcher_t::neighbours_of(index_type const ni) const
{
    meshless_sph_node_t const& node = (*nodes_)[ni];
    scalar_type const r             = node.kernel().h();
    Eigen::Vector3d const x         = node.kernel().xi();
    Discregrid::BoundingSphere const node_shape_function_domain{x, r};

    auto const intersects =
        [this, node_shape_function_domain](unsigned int node_idx, unsigned int depth) -> bool {
        Discregrid::BoundingSphere const& s = this->hull(node_idx);
        return s.overlaps(node_shape_function_domain);
    };

    std::vector<index_type> neighbours{};
    auto const get_neighbours = [this, node_shape_function_domain, &neighbours, ni](
                                    unsigned int node_idx,
                                    unsigned int depth) {
        base_type::Node const& node = this->node(node_idx);
        if (!node.isLeaf())
            return;

        for (auto j = node.begin; j < node.begin + node.n; ++j)
        {
            index_type const nj = m_lst[j];
            if (nj == ni)
                continue;

            Eigen::Vector3d const& xj = (*nodes_)[nj].Xi();
            if (node_shape_function_domain.contains(xj))
            {
                neighbours.push_back(nj);
            }
        }
    };

    traverseBreadthFirst(intersects, get_neighbours);
    neighbours.push_back(ni);
    return neighbours;
}

std::vector<index_type> meshless_sph_body_range_searcher_t::neighbours_of(
    Eigen::Vector3d const& p,
    scalar_type const h) const
{
    Discregrid::BoundingSphere const range{p, h};
    auto const intersects = [this, range](unsigned int node_idx, unsigned int depth) -> bool {
        Discregrid::BoundingSphere const& s = this->hull(node_idx);
        return s.overlaps(range);
    };

    std::vector<index_type> neighbours{};
    auto const get_neighbours =
        [this, range, &neighbours](unsigned int node_idx, unsigned int depth) {
            base_type::Node const& node = this->node(node_idx);
            if (!node.isLeaf())
                return;

            for (auto j = node.begin; j < node.begin + node.n; ++j)
            {
                index_type const nj = m_lst[j];

                Eigen::Vector3d const& xj = (*nodes_)[nj].Xi();
                if (range.contains(xj))
                {
                    neighbours.push_back(nj);
                }
            }
        };

    traverseBreadthFirst(intersects, get_neighbours);
    return neighbours;
}

std::vector<index_type>
meshless_sph_body_range_searcher_t::is_in_node_domains(Eigen::Vector3d const& p) const
{
    auto const intersects = [this, p](unsigned int node_idx, unsigned int depth) -> bool {
        Discregrid::BoundingSphere const& s = this->hull(node_idx);
        return s.contains(p);
    };

    std::vector<index_type> meshless_nodes{};
    auto const get_nodes = [this, p, &meshless_nodes](unsigned int node_idx, unsigned int depth) {
        base_type::Node const& node = this->node(node_idx);
        if (!node.isLeaf())
            return;

        for (auto j = node.begin; j < node.begin + node.n; ++j)
        {
            index_type const nj = m_lst[j];

            Eigen::Vector3d const& xj = (*nodes_)[nj].Xi();
            scalar_type const h       = (*nodes_)[nj].kernel().h();
            Discregrid::BoundingSphere const s(xj, h);
            if (s.contains(p))
            {
                meshless_nodes.push_back(nj);
            }
        }
    };

    traverseBreadthFirst(intersects, get_nodes);
    return meshless_nodes;
}

Eigen::Vector3d meshless_sph_body_range_searcher_t::entityPosition(unsigned int i) const
{
    return (*nodes_)[i].Xi();
}

void meshless_sph_body_range_searcher_t::computeHull(
    unsigned int b,
    unsigned int n,
    Discregrid::BoundingSphere& hull) const
{
    auto vertices_of_sphere = std::vector<Eigen::Vector3d>(n);
    for (unsigned int i = b; i < n + b; ++i)
        vertices_of_sphere[i - b] = (*nodes_)[m_lst[i]].Xi();

    Discregrid::BoundingSphere const s(vertices_of_sphere);

    hull.x() = s.x();
    hull.r() = s.r();
}

} // namespace mechanics
} // namespace physics
} // namespace sbs