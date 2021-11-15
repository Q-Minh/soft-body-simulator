#ifndef SBS_PHYSICS_VISUAL_MESHLESS_EMBEDDED_SURFACE_H
#define SBS_PHYSICS_VISUAL_MESHLESS_EMBEDDED_SURFACE_H

#include "interpolated_embedded_surface.h"

namespace sbs {
namespace physics {
namespace visual {

template <class MeshlessModelType>
class meshless_embedded_surface_t : public interpolated_embedded_surface_t<
                                        typename MeshlessModelType::interpolation_function_type>
{
  public:
    using meshless_model_type         = MeshlessModelType;
    using interpolation_function_type = typename MeshlessModelType::interpolation_function_type;
    using base_type = interpolated_embedded_surface_t<interpolation_function_type>;

    meshless_embedded_surface_t() = default;
    meshless_embedded_surface_t(
        std::vector<Eigen::Vector3d> const& points,
        std::vector<index_type> const& indices,
        meshless_model_type const* mechanical_model);

    meshless_model_type const* mechanical_model() const;
    std::vector<index_type> const& meshless_neighbours_of_vertex(index_type const vi) const;

    /**
     * @brief
     * Computes new vertex positions based on the mechanical model's current interpolation field
     */
    void update();

  private:
    meshless_model_type const* mechanical_model_;
    std::vector<std::vector<index_type>>
        neighbours_of_; ///< Precomputed neighbourhoods of each surface vertex to each surrounding
                        ///< meshless node
};

template <class MeshlessModelType>
inline meshless_embedded_surface_t<MeshlessModelType>::meshless_embedded_surface_t(
    std::vector<Eigen::Vector3d> const& points,
    std::vector<index_type> const& indices,
    meshless_model_type const* mechanical_model)
    : base_type(points, indices), mechanical_model_(mechanical_model)
{
    std::vector<interpolation_function_type> interpolation_functions{};
    interpolation_functions.reserve(this->vertex_count());

    auto const& Xs = this->reference_positions();
    for (auto vi = 0u; vi < Xs.size(); ++vi)
    {
        auto const& X = Xs[vi];
        interpolation_function_type const interpolation_function =
            this->mechanical_model_->interpolation_field_at(X);
        interpolation_functions.push_back(interpolation_function);

        // Get neighbours of surface vertex vi
        std::vector<index_type> const neighbours = this->mechanical_model_->in_support_of_nodes(X);
        neighbours_of_.push_back(neighbours);
    }

    this->use_interpolation_operators(interpolation_functions);
    for (auto i = 0u; i < Xs.size(); ++i)
    {
        auto const& Xi    = Xs[i];
        this->position(i) = Xi;
    }
}

template <class MeshlessModelType>
inline MeshlessModelType const*
meshless_embedded_surface_t<MeshlessModelType>::mechanical_model() const
{
    return mechanical_model_;
}

template <class MeshlessModelType>
inline std::vector<index_type> const&
meshless_embedded_surface_t<MeshlessModelType>::meshless_neighbours_of_vertex(
    index_type const vi) const
{
    return this->neighbours_of_[vi];
}

template <class MeshlessModelType>
inline void meshless_embedded_surface_t<MeshlessModelType>::update()
{
    auto const& Xs = this->reference_positions();
    for (auto i = 0u; i < Xs.size(); ++i)
    {
        auto const& Xi    = Xs[i];
        auto& interpolate = this->interpolation_operator(i);
        auto const xi     = interpolate(Xi);
        this->position(i) = xi.cast<scalar_type>();
    }
}

} // namespace visual
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_VISUAL_MESHLESS_EMBEDDED_SURFACE_H