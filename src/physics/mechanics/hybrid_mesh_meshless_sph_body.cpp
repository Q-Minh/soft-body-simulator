#include <Discregrid/cubic_lagrange_discrete_grid.hpp>
#include <Discregrid/geometry/mesh_distance.hpp>
#include <Discregrid/mesh/triangle_mesh.hpp>
#include <cassert>
#include <sbs/common/geometry.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_body.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_node.h>
#include <sbs/physics/particle.h>
#include <sbs/physics/simulation.h>

namespace sbs {
namespace physics {
namespace mechanics {

hybrid_mesh_meshless_mls_body_t::hybrid_mesh_meshless_mls_body_t(
    simulation_t& simulation,
    index_type id,
    common::geometry_t const& geometry,
    scalar_type const h,
    std::array<unsigned int, 3u> const& resolution)
    : body_t(simulation, id),
      meshless_nodes_(),
      material_space_meshless_node_searcher_(),
      volumetric_topology_(),
      is_boundary_tetrahedron_(),
      is_boundary_vertex_(),
      Ainv_(),
      visual_model_(),
      collision_model_(),
      h_(),
      mesh_particles_index_offset_(),
      meshless_particles_index_offset_()
{
    assert(geometry.geometry_type == common::geometry_t::geometry_type_t::tetrahedron);
    assert(geometry.has_indices());
    assert(geometry.has_positions());

    /**
     * Build tetrahedral exterior_layer_mesh topology
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
     * Create mesh meshless_nodes for fem shape functions. These mesh shape functions
     * will be precomputed in the initialize_physical_model() method
     */
    mesh_x0_.reserve(volumetric_topology_.vertex_count());
    is_boundary_vertex_.reserve(volumetric_topology_.vertex_count());
    for (std::size_t i = 0u; i < volumetric_topology_.vertex_count(); ++i)
    {
        auto const idx      = i * 3u;
        scalar_type const x = static_cast<scalar_type>(geometry.positions[idx]);
        scalar_type const y = static_cast<scalar_type>(geometry.positions[idx + 1u]);
        scalar_type const z = static_cast<scalar_type>(geometry.positions[idx + 2u]);
        Eigen::Vector3d const pos{x, y, z};
        mesh_x0_.push_back(pos);
        mesh_x_.push_back(pos);

        index_type const vi = static_cast<index_type>(i);
        // boundary tet mesh mesh_x0_ do not have associated shape functions,
        // so we do not create any XPBD particle for them
        bool const is_boundary_vertex = volumetric_topology_.is_boundary_vertex(vi);
        is_boundary_vertex_.push_back(is_boundary_vertex);
    }

    /**
     * Extract boundary surface of the tetrahedral exterior_layer_mesh
     */
    std::vector<std::array<unsigned int, 3u>> exterior_boundary_faces{};
    std::vector<triangle_t> const exterior_layer_boundary_triangles =
        volumetric_topology_.oriented_boundary_triangles();
    exterior_boundary_faces.reserve(exterior_layer_boundary_triangles.size());
    for (triangle_t const& f : exterior_layer_boundary_triangles)
        exterior_boundary_faces.push_back({f.v1(), f.v2(), f.v3()});

    /**
     * Extract interior layer's boundary surface of the tetrahedral exterior_layer_mesh
     */
    tetrahedron_set_t interior_tetrahedral_topology = volumetric_topology_;
    is_boundary_tetrahedron_.resize(volumetric_topology_.tetrahedron_count(), false);
    std::vector<index_type> const boundary_tets =
        volumetric_topology_.boundary_tetrahedron_indices();
    for (index_type const ti : boundary_tets)
    {
        is_boundary_tetrahedron_[ti] = true;
        interior_tetrahedral_topology.remove_tetrahedron(ti);
    }
    interior_tetrahedral_topology.collect_garbage();

    std::vector<std::array<unsigned int, 3u>> interior_boundary_faces{};
    std::vector<triangle_t> const interior_layer_boundary_triangles =
        interior_tetrahedral_topology.oriented_boundary_triangles();
    interior_boundary_faces.reserve(interior_layer_boundary_triangles.size());
    for (triangle_t const& f : interior_layer_boundary_triangles)
        interior_boundary_faces.push_back({f.v1(), f.v2(), f.v3()});

    /**
     * Compute the SDF of the exterior boundary surface and of the interior boundary surface to
     * obtain an inside/outside predicate for our particle sampling
     */
    Discregrid::TriangleMesh exterior_layer_mesh(mesh_x0_, exterior_boundary_faces);
    Discregrid::TriangleMesh interior_layer_mesh(mesh_x0_, interior_boundary_faces);
    Eigen::AlignedBox3d domain{};
    for (auto const& x : exterior_layer_mesh.vertices())
    {
        for (auto const& x : exterior_layer_mesh.vertices())
        {
            domain.extend(x);
        }
    }

    // add a bit of tolerance to the domain to fully enclose the exterior_layer_mesh
    domain.min() -= Eigen::Vector3d{1e-3, 1e-3, 1e-3};
    domain.max() += Eigen::Vector3d{1e-3, 1e-3, 1e-3};

    Discregrid::MeshDistance exterior_layer_md(exterior_layer_mesh);
    Discregrid::MeshDistance interior_layer_md(interior_layer_mesh);

    auto const exterior_layer_sdf = [&exterior_layer_md](Eigen::Vector3d const& xi) {
        return exterior_layer_md.signedDistanceCached(xi);
    };
    auto const interior_layer_sdf = [&interior_layer_md](Eigen::Vector3d const& xi) {
        return interior_layer_md.signedDistanceCached(xi);
    };

    Discregrid::CubicLagrangeDiscreteGrid grid(domain, resolution);
    auto const exterior_sdf_id = grid.addFunction(exterior_layer_sdf);
    auto const interior_sdf_id = grid.addFunction(interior_layer_sdf);

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
    // that lies inside the exterior boundary surface and outside the interior
    // boundary surface using the SDFs.
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

                scalar_type const exterior_layer_signed_distance =
                    grid.interpolate(exterior_sdf_id, grid_node_position);
                scalar_type const interior_layer_signed_distance =
                    grid.interpolate(interior_sdf_id, grid_node_position);

                bool const is_grid_node_inside_exterior_boundary =
                    exterior_layer_signed_distance < 0.;
                bool const is_grid_node_outside_interior_boundary =
                    interior_layer_signed_distance > 0.;
                bool const should_sample_particle =
                    is_grid_node_inside_exterior_boundary && is_grid_node_outside_interior_boundary;

                if (!should_sample_particle)
                    continue;

                functions::poly6_kernel_t kernel(grid_node_position, h_);
                auto const ni = static_cast<index_type>(meshless_nodes_.size());
                hybrid_mesh_meshless_mls_node_t node(ni, kernel, *this);
                node.xi() = kernel.xi();
                meshless_nodes_.push_back(node);
            }
        }
    }

    std::vector<Eigen::Vector3d> surface_vertices{};
    std::vector<triangle_t> surface_triangles{};
    std::unordered_map<index_type, index_type> tet_to_surface_vertex_map{};
    for (auto const& f : exterior_layer_boundary_triangles)
    {
        if (tet_to_surface_vertex_map.find(f.v1()) == tet_to_surface_vertex_map.end())
        {
            auto const index                  = tet_to_surface_vertex_map.size();
            tet_to_surface_vertex_map[f.v1()] = static_cast<index_type>(index);
            surface_vertices.push_back(mesh_x0_[f.v1()]);
        }
        if (tet_to_surface_vertex_map.find(f.v2()) == tet_to_surface_vertex_map.end())
        {
            auto const index                  = tet_to_surface_vertex_map.size();
            tet_to_surface_vertex_map[f.v2()] = static_cast<index_type>(index);
            surface_vertices.push_back(mesh_x0_[f.v2()]);
        }
        if (tet_to_surface_vertex_map.find(f.v3()) == tet_to_surface_vertex_map.end())
        {
            auto const index                  = tet_to_surface_vertex_map.size();
            tet_to_surface_vertex_map[f.v3()] = static_cast<index_type>(index);
            surface_vertices.push_back(mesh_x0_[f.v3()]);
        }

        auto const v1 = tet_to_surface_vertex_map[f.v1()];
        auto const v2 = tet_to_surface_vertex_map[f.v2()];
        auto const v3 = tet_to_surface_vertex_map[f.v3()];
        triangle_t const triangle{v1, v2, v3};
        surface_triangles.push_back(triangle);
    }

    // Compute the initial visual mesh as the boundary surface mesh of the given geometry
    visual_model_ = hybrid_mesh_meshless_sph_surface_t(this, surface_vertices, surface_triangles);
    for (std::size_t i = 0u; i < mesh_x0_.size(); ++i)
    {
        if (tet_to_surface_vertex_map.find(static_cast<index_type>(i)) ==
            tet_to_surface_vertex_map.end())
            continue;

        auto const idx = i * 3u;
        float const r  = static_cast<float>(geometry.colors[idx] / 255.f);
        float const g  = static_cast<float>(geometry.colors[idx + 1u] / 255.f);
        float const b  = static_cast<float>(geometry.colors[idx + 2u] / 255.f);

        auto const vi = tet_to_surface_vertex_map[static_cast<index_type>(i)];
        visual_model_.world_space_vertex(vi).color = Eigen::Vector3f{r, g, b};
    }
}

body_t::visual_model_type const& hybrid_mesh_meshless_mls_body_t::visual_model() const
{
    return visual_model_;
}

body_t::collision_model_type const& hybrid_mesh_meshless_mls_body_t::collision_model() const
{
    return collision_model_;
}

body_t::visual_model_type& hybrid_mesh_meshless_mls_body_t::visual_model()
{
    return visual_model_;
}

body_t::collision_model_type& hybrid_mesh_meshless_mls_body_t::collision_model()
{
    return collision_model_;
}

void hybrid_mesh_meshless_mls_body_t::update_visual_model()
{
    visual_model_.compute_positions();
    visual_model_.compute_normals();
}

void hybrid_mesh_meshless_mls_body_t::update_collision_model()
{
    collision_model_.update(simulation());
}

void hybrid_mesh_meshless_mls_body_t::update_physical_model()
{
    auto const& particles = simulation().particles()[id()];
    assert(particles.size() == (mesh_x0_.size() + meshless_nodes_.size()));
    assert(mesh_x_.size() == mesh_x0_.size());
    for (std::size_t i = 0u; i < mesh_x_.size(); ++i)
    {
        particle_t const& p = particles[mesh_particles_index_offset_ + i];
        mesh_x_[i]          = p.x();
    }
    for (std::size_t i = 0u; i < meshless_nodes_.size(); ++i)
    {
        index_type const ni                   = static_cast<index_type>(i);
        hybrid_mesh_meshless_mls_node_t& node = meshless_nodes_[ni];
        particle_t const& p                   = particles[meshless_particles_index_offset_ + i];
        node.xi()                             = p.x();
    }
}

void hybrid_mesh_meshless_mls_body_t::transform(Eigen::Affine3d const& affine)
{
    auto const v             = Eigen::Vector3d{1., 1., 1.}.normalized();
    auto const vp            = affine * Eigen::Vector4d{v.x(), v.y(), v.z(), 0.};
    auto const length_change = vp.norm();
    h_ *= length_change;

    for (std::size_t i = 0u; i < mesh_x_.size(); ++i)
    {
        mesh_x_[i]  = affine * mesh_x_[i].homogeneous();
        mesh_x0_[i] = affine * mesh_x0_[i].homogeneous();
    }
    for (std::size_t i = 0u; i < visual_model_.vertex_count(); ++i)
    {
        visual_model_.material_space_position(i) =
            affine * visual_model_.material_space_position(i).homogeneous();
        visual_model_.world_space_vertex(i).position =
            affine * visual_model_.world_space_vertex(i).position.homogeneous();
    }
    for (auto& meshless_node : meshless_nodes_)
    {
        meshless_node.Xi()         = affine * meshless_node.Xi().homogeneous();
        meshless_node.xi()         = affine * meshless_node.xi().homogeneous();
        meshless_node.kernel().h() = h_;
    }
}

void hybrid_mesh_meshless_mls_body_t::initialize_physical_model()
{
    // Precompute the mesh node shape functions as a function which maps
    // material space positions to barycentric coordinates. The barycentric
    // coordinates correspond to the shape functions of each tetrahedron vertex.
    Ainv_.reserve(volumetric_topology_.tetrahedron_count());
    for (std::size_t i = 0u; i < volumetric_topology_.tetrahedron_count(); ++i)
    {
        index_type const ti       = static_cast<index_type>(i);
        tetrahedron_t const& t    = volumetric_topology_.tetrahedron(ti);
        Eigen::Vector3d const& X0 = mesh_x0_[t.v1()];
        Eigen::Vector3d const& X1 = mesh_x0_[t.v2()];
        Eigen::Vector3d const& X2 = mesh_x0_[t.v3()];
        Eigen::Vector3d const& X3 = mesh_x0_[t.v4()];

        Eigen::Matrix4d A{};
        A.row(0u).setOnes();
        A.block(1, 0, 3, 1) = X0;
        A.block(1, 1, 3, 1) = X1;
        A.block(1, 2, 3, 1) = X2;
        A.block(1, 3, 3, 1) = X3;

        bool invertible{false};
        Eigen::Matrix4d Ainv{};
        double constexpr error = 1e-18;
        A.computeInverseWithCheck(Ainv, invertible, error);
        assert(invertible);
        Ainv_.push_back(Ainv);
    }

    // Get a spatial acceleration query object to obtain neighbours of
    // each meshless node in material space
    material_space_meshless_node_searcher_ =
        detail::hybrid_mesh_meshless_sph::meshless_node_range_searcher_t(&meshless_nodes_);

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
            material_space_meshless_node_searcher_.neighbours_of(ni);

        // precompute all neighbour information required to initialize our meshless meshless_nodes
        for (std::size_t k = 0u; k < neighbours_in_domain.size(); ++k)
        {
            index_type const j                   = neighbours_in_domain[k];
            meshless_sph_node_t const& neighbour = meshless_nodes_[j];
            // neighbour positions
            Eigen::Vector3d const* Xj = &neighbour.Xi();
            Xjs[i].push_back(Xj);
            // neighbour indices
            node_neighbour_indices[i].push_back(j);
            // neighbour meshless meshless_nodes
            node_neighbours[i].push_back(&neighbour);
        }
    }

    material_space_mesh_node_searcher_ =
        detail::hybrid_mesh_meshless_sph::mesh_tetrahedron_range_searcher_t(
            &volumetric_topology_,
            &mesh_x0_,
            &Ainv_);

    // update meshless meshless_nodes using precomputed neighbour information
    for (std::size_t i = 0u; i < meshless_nodes_.size(); ++i)
    {
        auto const& Xi      = meshless_nodes_[i].Xi();
        index_type const ti = material_space_mesh_node_searcher_.in_tetrahedron(Xi);

        // Used for debugging purposes only. Since we are testing inside/outside against
        // the cubic lagrange interpolated SDF at the SDF grid's resolution, some points
        // might be evaluated as "inside" the SDF, but are not actually exactly inside the
        // mesh.
        bool const has_found_parent = ti != std::numeric_limits<index_type>::max();

        auto const& xi = meshless_nodes_[i].xi();
        meshless_nodes_[i].initialize(xi, Xjs[i], node_neighbour_indices[i], ti);
    }
    // precompute the correction matrix Li for each meshless node
    for (std::size_t i = 0u; i < meshless_nodes_.size(); ++i)
    {
        meshless_nodes_[i].cache_Li_Vj(node_neighbours[i]);
    }
}

void hybrid_mesh_meshless_mls_body_t::initialize_visual_model()
{
    visual_model_.initialize_interpolation_scheme(2. * h_);
    visual_model_.compute_normals();
}

void hybrid_mesh_meshless_mls_body_t::initialize_collision_model()
{
    collision_model_      = collision::point_bvh_model_t(&visual_model_);
    collision_model_.id() = id();
}

void hybrid_mesh_meshless_mls_body_t::set_mesh_particles_index_offset(index_type offset)
{
    mesh_particles_index_offset_ = offset;
}

void hybrid_mesh_meshless_mls_body_t::set_meshless_particles_index_offset(index_type offset)
{
    meshless_particles_index_offset_ = offset;
}

index_type hybrid_mesh_meshless_mls_body_t::get_mesh_particles_index_offset() const
{
    return mesh_particles_index_offset_;
}

index_type hybrid_mesh_meshless_mls_body_t::get_meshless_particles_index_offset() const
{
    return meshless_particles_index_offset_;
}

std::vector<hybrid_mesh_meshless_mls_node_t> const&
hybrid_mesh_meshless_mls_body_t::meshless_nodes() const
{
    return meshless_nodes_;
}

std::vector<hybrid_mesh_meshless_mls_node_t>& hybrid_mesh_meshless_mls_body_t::meshless_nodes()
{
    return meshless_nodes_;
}

std::size_t hybrid_mesh_meshless_mls_body_t::mesh_node_count() const
{
    return mesh_x0_.size();
}

std::size_t hybrid_mesh_meshless_mls_body_t::meshless_node_count() const
{
    return meshless_nodes_.size();
}

hybrid_mesh_meshless_sph_surface_t const& hybrid_mesh_meshless_mls_body_t::surface_mesh() const
{
    return visual_model_;
}

hybrid_mesh_meshless_sph_surface_t& hybrid_mesh_meshless_mls_body_t::surface_mesh()
{
    return visual_model_;
}

collision::point_bvh_model_t const& hybrid_mesh_meshless_mls_body_t::bvh() const
{
    return collision_model_;
}

std::size_t hybrid_mesh_meshless_mls_body_t::mixed_meshless_node_count() const
{
    auto const count = std::count_if(
        meshless_nodes_.begin(),
        meshless_nodes_.end(),
        [](hybrid_mesh_meshless_mls_node_t const& meshless_node) {
            return meshless_node.is_mixed_particle();
        });
    return count;
}

std::size_t hybrid_mesh_meshless_mls_body_t::interior_tetrahedron_count() const
{
    auto const count =
        std::count(is_boundary_tetrahedron_.begin(), is_boundary_tetrahedron_.end(), false);
    return count;
}

std::size_t hybrid_mesh_meshless_mls_body_t::mesh_shape_function_count() const
{
    auto const count = std::count(is_boundary_vertex_.begin(), is_boundary_vertex_.end(), false);
    return count;
}

tetrahedron_set_t const& hybrid_mesh_meshless_mls_body_t::topology() const
{
    return volumetric_topology_;
}

tetrahedron_set_t& hybrid_mesh_meshless_mls_body_t::topology()
{
    return volumetric_topology_;
}

scalar_type hybrid_mesh_meshless_mls_body_t::h() const
{
    return h_;
}

std::vector<Eigen::Vector3d> const& hybrid_mesh_meshless_mls_body_t::x0() const
{
    return mesh_x0_;
}

std::vector<Eigen::Vector3d> const& hybrid_mesh_meshless_mls_body_t::x() const
{
    return mesh_x_;
}

detail::hybrid_mesh_meshless_sph::mesh_tetrahedron_range_searcher_t const&
hybrid_mesh_meshless_mls_body_t::mesh_tetrahedron_range_searcher() const
{
    return material_space_mesh_node_searcher_;
}

detail::hybrid_mesh_meshless_sph::meshless_node_range_searcher_t const&
hybrid_mesh_meshless_mls_body_t::meshless_node_range_searcher() const
{
    return material_space_meshless_node_searcher_;
}

bool hybrid_mesh_meshless_mls_body_t::is_boundary_mesh_tetrahedron(index_type const ti) const
{
    return is_boundary_tetrahedron_[ti];
}

bool hybrid_mesh_meshless_mls_body_t::is_boundary_mesh_vertex(index_type const vi) const
{
    return is_boundary_vertex_[vi];
}

Eigen::Matrix4d const& hybrid_mesh_meshless_mls_body_t::phi_i(index_type const ti) const
{
    return Ainv_[ti];
}

Eigen::Vector4d hybrid_mesh_meshless_mls_body_t::phi_i(index_type const ti, std::uint8_t v) const
{
    return Ainv_[ti].row(v);
}

Eigen::Matrix<scalar_type, 4, 3>
hybrid_mesh_meshless_mls_body_t::grad_phi_i(index_type const ti) const
{
    Eigen::Matrix4d const& Ainv = Ainv_[ti];
    Eigen::Matrix<scalar_type, 4, 3> grad{};
    for (auto i = 0u; i < 4u; ++i)
    {
        grad.row(i) = Ainv.block(i, 1u, 1u, 3u);
    }
    return grad;
}

Eigen::Vector3d
hybrid_mesh_meshless_mls_body_t::grad_phi_i(index_type const ti, std::uint8_t v) const
{
    Eigen::Matrix4d const& Ainv = Ainv_[ti];
    Eigen::Vector3d const grad  = Ainv.block(v, 1u, 1u, 3u).transpose();
    return grad;
}

/**
 * Spatial acceleration structures for quickly searching meshless meshless_nodes as well as mesh
 * meshless_nodes
 */
namespace detail {
namespace hybrid_mesh_meshless_sph {

meshless_node_range_searcher_t::meshless_node_range_searcher_t() : base_type(0u), nodes_() {}

meshless_node_range_searcher_t::meshless_node_range_searcher_t(
    std::vector<hybrid_mesh_meshless_mls_node_t> const* nodes)
    : base_type(nodes->size()), nodes_(nodes)
{
    this->construct();
}

std::vector<index_type> meshless_node_range_searcher_t::neighbours_of(index_type const ni) const
{
    hybrid_mesh_meshless_mls_node_t const& node = (*nodes_)[ni];
    scalar_type const r                         = node.kernel().h();
    Eigen::Vector3d const x                     = node.kernel().xi();
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

std::vector<index_type>
meshless_node_range_searcher_t::neighbours_of(Eigen::Vector3d const& p, scalar_type const h) const
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

Eigen::Vector3d meshless_node_range_searcher_t::entityPosition(unsigned int i) const
{
    return (*nodes_)[i].Xi();
}

void meshless_node_range_searcher_t::computeHull(
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

mesh_tetrahedron_range_searcher_t::mesh_tetrahedron_range_searcher_t()
    : base_type(0u), topology_(), mesh_nodes_(), Ainv_()
{
}

mesh_tetrahedron_range_searcher_t::mesh_tetrahedron_range_searcher_t(
    tetrahedron_set_t const* topology,
    std::vector<Eigen::Vector3d> const* mesh_nodes,
    std::vector<Eigen::Matrix4d> const* Ainv)
    : base_type(topology->tetrahedron_count()),
      topology_(topology),
      mesh_nodes_(mesh_nodes),
      Ainv_(Ainv)
{
    tetrahedron_centers_.reserve(topology_->tetrahedron_count());
    for (auto const& t : topology_->tetrahedra())
    {
        Eigen::Vector3d const& p1 = (*mesh_nodes)[t.v1()];
        Eigen::Vector3d const& p2 = (*mesh_nodes)[t.v2()];
        Eigen::Vector3d const& p3 = (*mesh_nodes)[t.v3()];
        Eigen::Vector3d const& p4 = (*mesh_nodes)[t.v4()];

        Eigen::Vector3d const center = 0.25 * (p1 + p2 + p3 + p4);
        tetrahedron_centers_.push_back(center);
    }
    this->construct();
}

index_type mesh_tetrahedron_range_searcher_t::in_tetrahedron(Eigen::Vector3d const& p) const
{
    auto const intersects = [this, p](unsigned int node_idx, unsigned int depth) -> bool {
        Discregrid::BoundingSphere const& s = this->hull(node_idx);
        return s.contains(p);
    };

    auto const is_point_in_tetrahedron =
        [this](Eigen::Vector3d const& point, tetrahedron_t const& t) {
            auto const& face_copies = t.faces_copy();
            std::array<bool, 4u> is_inside{false, false, false, false};
            for (std::uint8_t i = 0u; i < 4u; ++i)
            {
                auto const& f             = face_copies[i];
                Eigen::Vector3d const& p1 = (*mesh_nodes_)[f.v1()];
                Eigen::Vector3d const& p2 = (*mesh_nodes_)[f.v2()];
                Eigen::Vector3d const& p3 = (*mesh_nodes_)[f.v3()];

                Eigen::Vector3d const n           = (p2 - p1).cross(p3 - p1).normalized();
                scalar_type const signed_distance = (point - p1).dot(n);
                bool const is_outside             = signed_distance > 0.;
                is_inside[i]                      = !is_outside;
            }
            bool const is_p_in_t = is_inside[0] && is_inside[1] && is_inside[2] && is_inside[3];
            return is_p_in_t;
        };

    index_type parent_ti = std::numeric_limits<index_type>::max();
    bool found{false};
    auto const get_ti = [this, p, &parent_ti, is_point_in_tetrahedron, &found](
                            unsigned int node_idx,
                            unsigned int depth) {
        if (found)
            return;

        base_type::Node const& node = this->node(node_idx);
        if (!node.isLeaf())
            return;

        for (auto j = node.begin; j < node.begin + node.n; ++j)
        {
            index_type const ti    = static_cast<index_type>(m_lst[j]);
            tetrahedron_t const& t = topology_->tetrahedron(ti);

            Eigen::Vector3d const& p1 = (*mesh_nodes_)[t.v1()];
            Eigen::Vector3d const& p2 = (*mesh_nodes_)[t.v2()];
            Eigen::Vector3d const& p3 = (*mesh_nodes_)[t.v3()];
            Eigen::Vector3d const& p4 = (*mesh_nodes_)[t.v4()];

            Eigen::Matrix4d const& Ainv = (*Ainv_)[ti];

            Eigen::Vector4d const alpha = Ainv * Eigen::Vector4d{1., p.x(), p.y(), p.z()};

            scalar_type const a1 = alpha(0u);
            scalar_type const a2 = alpha(1u);
            scalar_type const a3 = alpha(2u);
            scalar_type const a4 = alpha(3u);

            scalar_type const alpha_sum = alpha.array().sum();

            if (is_point_in_tetrahedron(p, t))
            {
                // IMPORTANT NOTE:
                // Use this variable for debugging purposes.
                // It is possible that a point lies exactly in 2 tetrahedra or more.
                // For example, a point lying on a shared face of 2 tetrahedra will
                // belong in both tets. In this implementation, for the moment,
                // we will choose to associate a point with the first tetrahedra
                // that we find.
                bool const is_point_shared_between_multiple_tetrahedra =
                    (parent_ti != std::numeric_limits<index_type>::max());
                parent_ti = ti;
                found     = true;
            }
        }
    };

    traverseBreadthFirst(intersects, get_ti);

    return parent_ti;
}

Eigen::Vector3d mesh_tetrahedron_range_searcher_t::entityPosition(unsigned int i) const
{
    return tetrahedron_centers_[i];
}

void mesh_tetrahedron_range_searcher_t::computeHull(
    unsigned int b,
    unsigned int n,
    Discregrid::BoundingSphere& hull) const
{
    std::vector<Eigen::Vector3d> vertices_of_sphere{};
    vertices_of_sphere.reserve(n * 4u);
    for (unsigned int i = b; i < n + b; ++i)
    {
        index_type const ti    = static_cast<index_type>(m_lst[i]);
        tetrahedron_t const& t = topology_->tetrahedron(ti);
        for (index_type const vi : t.vertex_indices())
        {
            Eigen::Vector3d const& pi = (*mesh_nodes_)[vi];
            vertices_of_sphere.push_back(pi);
        }
    }

    Discregrid::BoundingSphere const s(vertices_of_sphere);

    hull.x() = s.x();
    hull.r() = s.r();
}

} // namespace hybrid_mesh_meshless_sph
} // namespace detail

} // namespace mechanics
} // namespace physics
} // namespace sbs
