#include <algorithm>
#include <sbs/common/geometry.h>
#include <sbs/physics/mechanics/meshless_body.h>
#include <sbs/physics/mechanics/meshless_node.h>
#include <sbs/physics/particle.h>
#include <sbs/physics/simulation.h>

namespace sbs {
namespace physics {
namespace mechanics {

meshless_body_t::meshless_body_t(
    simulation_t& simulation,
    index_type id,
    common::geometry_t const& geometry,
    scalar_type const h)
    : body_t(simulation, id),
      physical_model_(),
      material_space_range_query_(),
      volumetric_topology_(),
      visual_model_(),
      collision_model_(),
      h_(h)
{
    assert(geometry.geometry_type == common::geometry_t::geometry_type_t::tetrahedron);
    assert(geometry.has_indices());
    assert(geometry.has_positions());

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
    visual_model_ = tetrahedral_mesh_boundary_t(&volumetric_topology_);

    for (std::size_t i = 0u; i < volumetric_topology_.vertex_count(); ++i)
    {
        auto const idx      = i * 3u;
        scalar_type const x = static_cast<scalar_type>(geometry.positions[idx]);
        scalar_type const y = static_cast<scalar_type>(geometry.positions[idx + 1u]);
        scalar_type const z = static_cast<scalar_type>(geometry.positions[idx + 2u]);
        Eigen::Vector3d const pos{x, y, z};

        particle_t const p{pos};
        simulation.add_particle(p, this->id());
    }

    for (std::size_t i = 0u; i < visual_model_.vertex_count(); ++i)
    {
        auto const particle_index = visual_model_.from_surface_vertex(i);

        auto const idx = i * 3u;
        float const r  = static_cast<float>(geometry.colors[idx] / 255.f);
        float const g  = static_cast<float>(geometry.colors[idx + 1u] / 255.f);
        float const b  = static_cast<float>(geometry.colors[idx + 2u] / 255.f);

        visual_model_.mutable_vertex(i).position =
            simulation.particles()[this->id()][particle_index].x();
        visual_model_.mutable_vertex(i).color = Eigen::Vector3f{r, g, b};
    }
    visual_model_.compute_normals();

    collision_model_      = collision::point_bvh_model_t(&visual_model_);
    collision_model_.id() = this->id();
}

body_t::visual_model_type const& meshless_body_t::visual_model() const
{
    return visual_model_;
}

body_t::collision_model_type const& meshless_body_t::collision_model() const
{
    return collision_model_;
}

body_t::visual_model_type& meshless_body_t::visual_model()
{
    return visual_model_;
}

body_t::collision_model_type& meshless_body_t::collision_model()
{
    return collision_model_;
}

void meshless_body_t::update_visual_model()
{
    auto const& particles = simulation().particles()[this->id()];
    for (std::size_t i = 0u; i < visual_model_.vertex_count(); ++i)
    {
        index_type const node_i                  = visual_model_.from_surface_vertex(i);
        visual_model_.mutable_vertex(i).position = particles[node_i].x();
    }
    visual_model_.compute_normals();
}

void meshless_body_t::update_collision_model()
{
    collision_model_.update(simulation());
}

void meshless_body_t::update_physical_model()
{
    // no-op
}

void meshless_body_t::transform(Eigen::Affine3d const& affine)
{
    auto& particles = simulation().particles().at(id());
    for (auto& p : particles)
    {
        p.x0() = affine * p.x0().homogeneous();
        p.xi() = affine * p.xi().homogeneous();
        p.xn() = affine * p.xn().homogeneous();
        p.x()  = affine * p.x().homogeneous();
    }

    update_visual_model();
    collision_model_      = collision::point_bvh_model_t(&visual_model_);
    collision_model_.id() = this->id();
}

std::vector<meshless_node_t> const& meshless_body_t::nodes() const
{
    return physical_model_;
}

physics::tetrahedral_mesh_boundary_t const& meshless_body_t::surface_mesh() const
{
    return visual_model_;
}

physics::tetrahedral_mesh_boundary_t& meshless_body_t::surface_mesh()
{
    return visual_model_;
}

collision::point_bvh_model_t const& meshless_body_t::bvh() const
{
    return collision_model_;
}

void meshless_body_t::initialize_physical_model(simulation_t const& simulation)
{
    physical_model_.clear();

    std::vector<particle_t> const& particles = simulation.particles()[this->id()];
    physical_model_.reserve(particles.size());
    for (std::size_t i = 0u; i < particles.size(); ++i)
    {
        particle_t const& p = particles[i];
        functions::poly6_kernel_t kernel(p.x0(), h_);
        auto const ni = static_cast<index_type>(i);
        meshless_node_t node(ni, kernel);
        physical_model_.push_back(node);
    }
    material_space_range_query_ = range_searcher_t(&physical_model_);
    std::vector<std::vector<Eigen::Vector3d const*>> Xjs{};
    std::vector<std::vector<index_type>> node_neighbour_indices{};
    std::vector<std::vector<meshless_node_t const*>> node_neighbours{};
    std::vector<functions::poly6_kernel_t> node_kernels{};

    Xjs.resize(physical_model_.size());
    node_neighbour_indices.resize(physical_model_.size());
    node_neighbours.resize(physical_model_.size());
    node_kernels.resize(physical_model_.size());

    for (std::size_t i = 0u; i < physical_model_.size(); ++i)
    {
        meshless_node_t const& node = physical_model_[i];
        auto const ni               = static_cast<index_type>(i);
        std::vector<index_type> const neighbours_in_domain =
            material_space_range_query_.neighbours_of(ni);

        node_kernels[i] = node.kernel();

        // precompute all neighbour information required to initialize our meshless nodes
        for (std::size_t k = 0u; k < neighbours_in_domain.size(); ++k)
        {
            index_type const j               = neighbours_in_domain[k];
            meshless_node_t const& neighbour = physical_model_[j];
            // neighbour positions
            Eigen::Vector3d const* Xj = &neighbour.Xi();
            Xjs[i].push_back(Xj);
            // neighbour indices
            node_neighbour_indices[i].push_back(j);
            // neighbour meshless nodes
            node_neighbours[i].push_back(&neighbour);
        }
    }
    for (std::size_t i = 0u; i < physical_model_.size(); ++i)
    {
        auto const ni = static_cast<index_type>(i);
        physical_model_[i] =
            meshless_node_t(ni, Xjs[i], node_neighbour_indices[i], node_kernels[i]);
    }
    for (std::size_t i = 0u; i < physical_model_.size(); ++i)
    {
        physical_model_[i].cache_Li_Vj(node_neighbours[i]);
    }
}

range_searcher_t::range_searcher_t() : base_type(0u), nodes_() {}

range_searcher_t::range_searcher_t(std::vector<meshless_node_t> const* nodes)
    : base_type(nodes->size()), nodes_(nodes)
{
    this->construct();
}

std::vector<index_type> range_searcher_t::neighbours_of(index_type const ni) const
{
    meshless_node_t const& node = (*nodes_)[ni];
    scalar_type const r         = node.kernel().h();
    Eigen::Vector3d const x     = node.kernel().xi();
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

Eigen::Vector3d range_searcher_t::entityPosition(unsigned int i) const
{
    return (*nodes_)[i].Xi();
}

void range_searcher_t::computeHull(unsigned int b, unsigned int n, Discregrid::BoundingSphere& hull)
    const
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