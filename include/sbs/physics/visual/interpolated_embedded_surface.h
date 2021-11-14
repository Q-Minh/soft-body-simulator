#ifndef SBS_PHYSICS_VISUAL_INTERPOLATED_EMBEDDED_SURFACE_H
#define SBS_PHYSICS_VISUAL_INTERPOLATED_EMBEDDED_SURFACE_H

#include "sbs/aliases.h"
#include "sbs/common/node.h"

#include <Eigen/Core>
#include <array>
#include <vector>

namespace sbs {
namespace physics {
namespace visual {

template <class InterpolationFunctionType>
class interpolated_embedded_surface_t : public common::renderable_node_t
{
  public:
    using interpolation_op_type = InterpolationFunctionType;

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
    std::vector<Eigen::Vector3d> const& positions() const { return x_; }
    std::vector<index_type> const& indices() const { return indices_; }

    interpolation_op_type const& interpolation_operator(index_type vi) const { return maps_[vi]; }
    Eigen::Vector3d const& reference_position(index_type vi) const { return X_[vi]; }
    Eigen::Vector3d const& position(index_type vi) const { return x_[vi]; }
    std::array<index_type, 3u> triangle(index_type fi) const;

    std::size_t vertex_count() const { return X_.size(); }
    std::size_t triangle_count() const { return indices_.size() / 3u; }

    // Mutators
    void
    use_interpolation_operators(std::vector<interpolation_op_type> const& interpolation_operators);

    virtual void prepare_vertices_for_rendering() override;
    virtual void prepare_indices_for_rendering() override;

  protected:
    std::vector<interpolation_op_type>& interpolation_operators() { maps_; }
    std::vector<Eigen::Vector3d>& reference_positions() { return X_; }
    std::vector<Eigen::Vector3d>& positions() { return x_; }
    std::vector<index_type>& indices() { return indices_; }

    interpolation_op_type& interpolation_operator(index_type vi) { return maps_[vi]; }
    Eigen::Vector3d& reference_position(index_type vi) { return X_[vi]; }
    Eigen::Vector3d& position(index_type vi) { return x_[vi]; }

  private:
    std::vector<interpolation_op_type>
        maps_; ///< Mapping from material space positions to world space positions
    std::vector<Eigen::Vector3d> X_;  ///< Material space positions of the surface mesh vertices
    std::vector<Eigen::Vector3d> x_;  ///< World space positions of the surface mesh vertices
    std::vector<index_type> indices_; ///< Surface mesh triangle vertex indices
};

template <class InterpolationFunctionType>
interpolated_embedded_surface_t<InterpolationFunctionType>::interpolated_embedded_surface_t(
    std::vector<interpolation_op_type> const& interpolation_functions,
    std::vector<Eigen::Vector3d> const& vertex_positions,
    std::vector<index_type> const& indices)
    : maps_(interpolation_functions), X_(vertex_positions), x_(), indices_(indices)
{
    x_.reserve(X_.size());
    for (auto i = 0u; i < X_.size(); ++i)
    {
        auto const& Xi           = X_[i];
        auto const& interpolate  = maps_[i];
        Eigen::Vector3d const xi = interpolate(Xi).cast<scalar_type>();
        x_.push_back(xi);
    }
}
template <class InterpolationFunctionType>
inline interpolated_embedded_surface_t<InterpolationFunctionType>::interpolated_embedded_surface_t(
    std::vector<Eigen::Vector3d> const& vertex_positions,
    std::vector<index_type> const& indices)
    : maps_(), X_(vertex_positions), x_(vertex_positions), indices_(indices)
{
}

template <class InterpolationFunctionType>
std::array<index_type, 3u>
interpolated_embedded_surface_t<InterpolationFunctionType>::triangle(index_type fi) const
{
    auto const offset = fi * 3u;
    std::array<index_type, 3u> const f{
        indices_[offset + 0u],
        indices_[offset + 1u],
        indices_[offset + 2u]};

    return f;
}

template <class InterpolationFunctionType>
inline void interpolated_embedded_surface_t<InterpolationFunctionType>::use_interpolation_operators(
    std::vector<interpolation_op_type> const& interpolation_operators)
{
    maps_ = interpolation_operators;
}

template <class InterpolationFunctionType>
void interpolated_embedded_surface_t<InterpolationFunctionType>::prepare_vertices_for_rendering()
{
    std::vector<float> cpu_buffer{};
    cpu_buffer.reserve(x_.size());
    for (auto i = 0u; i < x_.size(); ++i)
    {
        cpu_buffer.push_back(static_cast<float>(x_[i].x()));
        cpu_buffer.push_back(static_cast<float>(x_[i].y()));
        cpu_buffer.push_back(static_cast<float>(x_[i].z()));
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
