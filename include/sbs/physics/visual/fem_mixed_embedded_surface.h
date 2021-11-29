#ifndef SBS_PHYSICS_VISUAL_FEM_MIXED_EMBEDDED_SURFACE_H
#define SBS_PHYSICS_VISUAL_FEM_MIXED_EMBEDDED_SURFACE_H

#include "interpolated_embedded_surface.h"

namespace sbs {
namespace physics {
namespace visual {

/**
 * @brief
 * Surface mesh embedded by interpolation in a mixed fem + meshless model.
 * For now, it is expected that the surface mesh resides entirely in the
 * mixed region (where interpolation is between fem + meshless basis functions).
 * @tparam FemMixedModelType
 */
template <class FemMixedModelType>
class fem_mixed_embedded_surface_t
    : public interpolated_embedded_surface_t<
          typename FemMixedModelType::mixed_interpolation_function_type>
{
  public:
    using fem_mixed_model_type = FemMixedModelType;
    using interpolation_function_type =
        typename fem_mixed_model_type::mixed_interpolation_function_type;
    using base_type = interpolated_embedded_surface_t<interpolation_function_type>;

    fem_mixed_embedded_surface_t() = default;
    fem_mixed_embedded_surface_t(
        std::vector<Eigen::Vector3d> const& points,
        std::vector<index_type> const& indices,
        fem_mixed_model_type* mechanical_model);

    fem_mixed_model_type const* mechanical_model() const;
    fem_mixed_model_type* mechanical_model();
    std::vector<index_type> const& neighbours_of(index_type const vi) const;
    index_type in_cell(index_type const vi) const;

    /**
     * @brief
     * Computes new vertex positions based on the mechanical model's current interpolation field
     */
    void update();

  private:
    fem_mixed_model_type* mechanical_model_;
    std::vector<index_type> in_cell_; ///< Englobing tetrahedra of surface vertices
    std::vector<std::vector<index_type>>
        neighbours_of_; ///< Precomputed neighbourhoods of each surface vertex to each surrounding
                        ///< meshless node
};

template <class FemMixedModelType>
inline fem_mixed_embedded_surface_t<FemMixedModelType>::fem_mixed_embedded_surface_t(
    std::vector<Eigen::Vector3d> const& points,
    std::vector<index_type> const& indices,
    fem_mixed_model_type* mechanical_model)
    : base_type(points, indices), mechanical_model_(mechanical_model), in_cell_(), neighbours_of_()
{
    std::vector<interpolation_function_type> interpolation_functions{};
    interpolation_functions.reserve(this->vertex_count());

    in_cell_.reserve(this->vertex_count());

    auto const& Xs = this->reference_positions();
    for (auto vi = 0u; vi < Xs.size(); ++vi)
    {
        auto const& X = Xs[vi];
        interpolation_function_type const interpolation_function =
            this->mechanical_model_->mixed_interpolation_field_at(X);
        interpolation_functions.push_back(interpolation_function);

        // Get englobing tetrahedron
        auto const& fem_domain = this->mechanical_model_->domain();
        index_type const e     = fem_domain.in_tetrahedron(X);
        in_cell_.push_back(e);

        // Get neighbours of surface vertex vi
        auto const& meshless_model               = this->mechanical_model_->meshless_model();
        std::vector<index_type> const neighbours = meshless_model.in_support_of_nodes(X);
        neighbours_of_.push_back(neighbours);
    }

    this->use_interpolation_operators(interpolation_functions);
    update();
}

template <class FemMixedModelType>
inline FemMixedModelType const*
fem_mixed_embedded_surface_t<FemMixedModelType>::mechanical_model() const
{
    return mechanical_model_;
}

template <class FemMixedModelType>
inline FemMixedModelType* fem_mixed_embedded_surface_t<FemMixedModelType>::mechanical_model()
{
    return mechanical_model_;
}

template <class FemMixedModelType>
inline std::vector<index_type> const&
fem_mixed_embedded_surface_t<FemMixedModelType>::neighbours_of(index_type const vi) const
{
    return this->neighbours_of_[vi];
}

template <class FemMixedModelType>
inline index_type
fem_mixed_embedded_surface_t<FemMixedModelType>::in_cell(index_type const vi) const
{
    return this->in_cell_[vi];
}

template <class FemMixedModelType>
inline void fem_mixed_embedded_surface_t<FemMixedModelType>::update()
{
    auto const& Xs = this->reference_positions();
    for (auto i = 0u; i < Xs.size(); ++i)
    {
        auto const& Xi          = Xs[i];
        auto const& interpolate = this->interpolation_operator(i);
        auto const xi           = interpolate(Xi);
        this->position(i)       = xi.cast<scalar_type>();
    }
    compute_normals();
}

} // namespace visual
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_VISUAL_FEM_MIXED_EMBEDDED_SURFACE_H
