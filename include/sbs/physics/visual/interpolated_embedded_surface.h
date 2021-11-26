#ifndef SBS_PHYSICS_VISUAL_INTERPOLATED_EMBEDDED_SURFACE_H
#define SBS_PHYSICS_VISUAL_INTERPOLATED_EMBEDDED_SURFACE_H

#include "sbs/aliases.h"
#include "sbs/common/mesh.h"

#include <Eigen/Core>
#include <array>
#include <vector>

namespace sbs {
namespace physics {
namespace visual {

template <class InterpolationFunctionType>
class interpolated_embedded_surface_t : public common::shared_vertex_surface_mesh_i
{
  public:
    using interpolation_op_type = InterpolationFunctionType;

    interpolated_embedded_surface_t() = default;

    /**
     * @brief
     * Creates an interpolated triangle surface mesh by
     * taking in the vertex positions of the mesh in some reference space,
     * and the topology as triangle indices. The world space vertex positions
     * of the triangle mesh can be recovered by interpolating each reference
     * space vertex position using its corresponding interpolation function.
     * @param interpolation_functions
     * @param vertex_positions
     * @param indices
     */
    interpolated_embedded_surface_t(
        std::vector<interpolation_op_type> const& interpolation_functions,
        std::vector<Eigen::Vector3d> const& vertex_positions,
        std::vector<index_type> const& indices);

    interpolated_embedded_surface_t(
        std::vector<Eigen::Vector3d> const& vertex_positions,
        std::vector<index_type> const& indices);

    // Accessors
    std::vector<interpolation_op_type> const& interpolation_operators() const { maps_; }
    std::vector<Eigen::Vector3d> const& reference_positions() const { return X_; }
    std::vector<vertex_type> const& vertices() const { return vertices_; }
    std::vector<index_type> const& indices() const { return indices_; }

    interpolation_op_type const& interpolation_operator(index_type vi) const { return maps_[vi]; }
    Eigen::Vector3d const& reference_position(index_type vi) const { return X_[vi]; }
    Eigen::Vector3d const& position(index_type vi) const { return vertices_[vi].position; }

    virtual std::size_t vertex_count() const override { return X_.size(); }
    virtual std::size_t triangle_count() const override { return indices_.size() / 3u; }

    virtual vertex_type vertex(std::size_t vi) const override { return vertices_[vi]; }
    virtual triangle_type triangle(std::size_t fi) const override;

    // Mutators
    void
    use_interpolation_operators(std::vector<interpolation_op_type> const& interpolation_operators);
    void set_color(Eigen::Vector3f const& rgb);
    void compute_normals();

    virtual void prepare_vertices_for_rendering() override;
    virtual void prepare_indices_for_rendering() override;

    std::vector<interpolation_op_type>& interpolation_operators() { maps_; }
    std::vector<Eigen::Vector3d>& reference_positions() { return X_; }
    std::vector<vertex_type>& vertices() { return vertices_; }
    std::vector<index_type>& indices() { return indices_; }

    interpolation_op_type& interpolation_operator(index_type vi) { return maps_[vi]; }
    Eigen::Vector3d& reference_position(index_type vi) { return X_[vi]; }
    Eigen::Vector3d& position(index_type vi) { return vertices_[vi].position; }

  private:
    std::vector<interpolation_op_type>
        maps_; ///< Mapping from material space positions to world space positions
    std::vector<Eigen::Vector3d> X_;    ///< Material space positions of the surface mesh vertices
    std::vector<vertex_type> vertices_; ///< World space vertices of the surface mesh
    std::vector<index_type> indices_;   ///< Surface mesh triangle vertex indices
};

template <class InterpolationFunctionType>
interpolated_embedded_surface_t<InterpolationFunctionType>::interpolated_embedded_surface_t(
    std::vector<interpolation_op_type> const& interpolation_functions,
    std::vector<Eigen::Vector3d> const& vertex_positions,
    std::vector<index_type> const& indices)
    : maps_(interpolation_functions), X_(vertex_positions), vertices_(), indices_(indices)
{
    vertices_.resize(X_.size());
    for (auto i = 0u; i < X_.size(); ++i)
    {
        auto const& Xi           = X_[i];
        auto const& interpolate  = maps_[i];
        Eigen::Vector3d const xi = interpolate(Xi).cast<scalar_type>();
        vertices_[i].position    = xi;
    }
    compute_normals();
}
template <class InterpolationFunctionType>
inline interpolated_embedded_surface_t<InterpolationFunctionType>::interpolated_embedded_surface_t(
    std::vector<Eigen::Vector3d> const& vertex_positions,
    std::vector<index_type> const& indices)
    : maps_(), X_(vertex_positions), vertices_(), indices_(indices)
{
    vertices_.resize(X_.size());
    for (auto i = 0u; i < X_.size(); ++i)
    {
        vertices_[i].position = X_[i];
    }
    compute_normals();
}

template <class InterpolationFunctionType>
inline common::shared_vertex_surface_mesh_i::triangle_type
interpolated_embedded_surface_t<InterpolationFunctionType>::triangle(std::size_t fi) const
{
    auto const offset = fi * 3u;
    std::array<index_type, 3u> const f{
        indices_[offset + 0u],
        indices_[offset + 1u],
        indices_[offset + 2u]};

    return triangle_type(f);
}

template <class InterpolationFunctionType>
inline void interpolated_embedded_surface_t<InterpolationFunctionType>::use_interpolation_operators(
    std::vector<interpolation_op_type> const& interpolation_operators)
{
    maps_ = interpolation_operators;
}

template <class InterpolationFunctionType>
inline void
interpolated_embedded_surface_t<InterpolationFunctionType>::set_color(Eigen::Vector3f const& rgb)
{
    for (auto& v : vertices_)
    {
        v.color = rgb;
    }
}

template <class InterpolationFunctionType>
inline void interpolated_embedded_surface_t<InterpolationFunctionType>::compute_normals()
{
    for (std::size_t i = 0u; i < indices_.size(); i += 3u)
    {
        auto const v1 = indices_[i + 0u];
        auto const v2 = indices_[i + 1u];
        auto const v3 = indices_[i + 2u];

        Eigen::Vector3d const& p1 = vertices_[v1].position;
        Eigen::Vector3d const& p2 = vertices_[v2].position;
        Eigen::Vector3d const& p3 = vertices_[v3].position;

        Eigen::Vector3d const n = (p2 - p1).cross(p3 - p1);

        vertices_[v1].normal += n;
        vertices_[v2].normal += n;
        vertices_[v3].normal += n;
    }

    for (std::size_t i = 0u; i < vertices_.size(); ++i)
    {
        vertices_[i].normal.normalize();
    }
}

template <class InterpolationFunctionType>
void interpolated_embedded_surface_t<InterpolationFunctionType>::prepare_vertices_for_rendering()
{
    std::vector<float> cpu_buffer{};
    cpu_buffer.reserve(vertices_.size() * 9u);
    for (auto const& v : vertices_)
    {
        cpu_buffer.push_back(static_cast<float>(v.position.x()));
        cpu_buffer.push_back(static_cast<float>(v.position.y()));
        cpu_buffer.push_back(static_cast<float>(v.position.z()));
        cpu_buffer.push_back(static_cast<float>(v.normal.x()));
        cpu_buffer.push_back(static_cast<float>(v.normal.y()));
        cpu_buffer.push_back(static_cast<float>(v.normal.z()));
        cpu_buffer.push_back(v.color.x());
        cpu_buffer.push_back(v.color.y());
        cpu_buffer.push_back(v.color.z());
    }
    this->transfer_vertices_for_rendering(std::move(cpu_buffer));
}

template <class InterpolationFunctionType>
void interpolated_embedded_surface_t<InterpolationFunctionType>::prepare_indices_for_rendering()
{
    std::vector<unsigned int> index_buffer{};
    index_buffer.reserve(indices_.size());
    std::copy(indices_.begin(), indices_.end(), std::back_inserter(index_buffer));
    this->transfer_indices_for_rendering(std::move(index_buffer));
}

} // namespace visual
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_VISUAL_INTERPOLATED_EMBEDDED_SURFACE_H
