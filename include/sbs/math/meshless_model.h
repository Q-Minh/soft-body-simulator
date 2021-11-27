#ifndef SBS_MATH_MESHLESS_MODEL_H
#define SBS_MATH_MESHLESS_MODEL_H

#include "sbs/aliases.h"

#include <Discregrid/acceleration/bounding_sphere.hpp>
#include <Discregrid/acceleration/kd_tree.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <vector>

namespace sbs {
namespace math {

template <class DofType, class PointType, class BasisFunctionType>
class meshless_model_t
{
  public:
    using self_type           = meshless_model_t<DofType, PointType, BasisFunctionType>;
    using dof_type            = DofType;
    using point_type          = PointType;
    using basis_function_type = BasisFunctionType;

    meshless_model_t() = default;
    meshless_model_t(meshless_model_t const& other);

    meshless_model_t& operator=(meshless_model_t const& other);

    // Accessors
    dof_type const& dof(index_type i) const { return dofs_[i]; }
    point_type const& point(index_type i) const { return points_[i]; }
    basis_function_type const& phi(index_type i) const { return phis_[i]; }

    dof_type& dof(index_type i) { return dofs_[i]; }
    basis_function_type& phi(index_type i) { return phis_[i]; }

    std::size_t dof_count() const { return dofs_.size(); }
    std::size_t point_count() const { return points_.size(); }

    scalar_type h() const { return h_; }
    std::vector<index_type> in_support_of_nodes(Eigen::Vector3d const& X) const;

    // Modifiers
    void add_dof(dof_type const& dof) { dofs_.push_back(dof); }
    void add_point(point_type const& point) { points_.push_back(point); }
    void add_basis_function(basis_function_type const& phi) { phis_.push_back(phi); }
    void set_support_radius(scalar_type h) { h_ = h; }

    void initialize_in_support_query(scalar_type tolerance);

    void clear_dofs() { dofs_.clear(); }
    void clear_points() { points_.clear(); }
    void clear()
    {
        clear_dofs();
        clear_points();
    }

    class in_support_query_t;

    std::vector<dof_type> const& dofs() const { return dofs_; }
    std::vector<dof_type>& dofs() { return dofs_; }

    std::vector<point_type> const& points() const { return points_; }
    std::vector<point_type>& points() { return points_; }

    std::vector<basis_function_type> const& phis() const { return phis_; }
    std::vector<basis_function_type>& phis() { return phis_; }

  private:
    std::vector<dof_type> dofs_;
    std::vector<point_type> points_;
    std::vector<basis_function_type> phis_;
    scalar_type h_;

  public:
    class in_support_query_t : public Discregrid::KDTree<Discregrid::BoundingSphere>
    {
      public:
        using base_type           = Discregrid::KDTree<Discregrid::BoundingSphere>;
        using meshless_model_type = meshless_model_t<dof_type, point_type, basis_function_type>;

        in_support_query_t();
        in_support_query_t(std::vector<point_type> const* points, scalar_type tolerance);

        in_support_query_t(in_support_query_t const& other) = default;
        in_support_query_t(in_support_query_t&& other)      = default;

        in_support_query_t& operator=(in_support_query_t const& other) = default;
        in_support_query_t& operator=(in_support_query_t&& other) noexcept = default;

        std::vector<index_type>
        neighbours_in_support(Eigen::Vector3d const& X, scalar_type radius) const;

        virtual Eigen::Vector3d entityPosition(unsigned int i) const override final;
        virtual void computeHull(unsigned int b, unsigned int n, Discregrid::BoundingSphere& hull)
            const override final;

        friend class meshless_model_t<DofType, PointType, BasisFunctionType>;

      private:
        std::vector<point_type> const* points_;
        scalar_type tolerance_;
    };

  private:
    in_support_query_t in_support_query_;
};

template <class DofType, class PointType, class BasisFunctionType>
inline meshless_model_t<DofType, PointType, BasisFunctionType>::meshless_model_t(
    meshless_model_t const& other)
    : dofs_(other.dofs_), points_(other.points_), phis_(other.phis_), h_(other.h_)
{
    in_support_query_ = in_support_query_t(&points_, other.in_support_query_.tolerance_);
}

template <class DofType, class PointType, class BasisFunctionType>
inline meshless_model_t<DofType, PointType, BasisFunctionType>&
meshless_model_t<DofType, PointType, BasisFunctionType>::operator=(meshless_model_t const& other)
{
    dofs_             = other.dofs_;
    points_           = other.points_;
    phis_             = other.phis_;
    h_                = other.h_;
    in_support_query_ = in_support_query_t(&points_, other.in_support_query_.tolerance_);
    return *this;
}

template <class DofType, class PointType, class BasisFunctionType>
inline std::vector<index_type>
meshless_model_t<DofType, PointType, BasisFunctionType>::in_support_of_nodes(
    Eigen::Vector3d const& X) const
{
    return in_support_query_.neighbours_in_support(X, h_);
}

template <class DofType, class PointType, class BasisFunctionType>
inline void meshless_model_t<DofType, PointType, BasisFunctionType>::initialize_in_support_query(
    scalar_type tolerance)
{
    in_support_query_ = in_support_query_t(&points_, tolerance);
}

template <class DofType, class PointType, class BasisFunctionType>
inline meshless_model_t<DofType, PointType, BasisFunctionType>::in_support_query_t::
    in_support_query_t()
    : base_type(0u), points_(nullptr), tolerance_(0.)
{
}

template <class DofType, class PointType, class BasisFunctionType>
inline meshless_model_t<DofType, PointType, BasisFunctionType>::in_support_query_t::
    in_support_query_t(std::vector<PointType> const* points, scalar_type tolerance)
    : base_type(points->size()), points_(points), tolerance_(tolerance)
{
    this->construct();
}

template <class DofType, class PointType, class BasisFunctionType>
inline std::vector<index_type>
meshless_model_t<DofType, PointType, BasisFunctionType>::in_support_query_t::neighbours_in_support(
    Eigen::Vector3d const& X,
    scalar_type radius) const
{
    Discregrid::BoundingSphere const support{X, radius};

    auto const intersects = [this, support](unsigned int node_idx, unsigned int depth) -> bool {
        Discregrid::BoundingSphere const& s = this->hull(node_idx);
        return s.overlaps(support);
    };

    std::vector<index_type> neighbours{};
    auto const get_neighbours =
        [this, support, &neighbours](unsigned int node_idx, unsigned int depth) {
            base_type::Node const& node = this->node(node_idx);
            if (!node.isLeaf())
                return;

            for (auto j = node.begin; j < node.begin + node.n; ++j)
            {
                index_type const nj = m_lst[j];

                Eigen::Vector3d const& Xj = (*points_)[nj].cast<scalar_type>();
                if (support.contains(Xj))
                {
                    neighbours.push_back(nj);
                }
            }
        };

    this->traverseBreadthFirst(intersects, get_neighbours);
    std::sort(neighbours.begin(), neighbours.end());
    auto it = std::unique(neighbours.begin(), neighbours.end());
    neighbours.erase(it, neighbours.end());
    return neighbours;
}

template <class DofType, class PointType, class BasisFunctionType>
inline Eigen::Vector3d
meshless_model_t<DofType, PointType, BasisFunctionType>::in_support_query_t::entityPosition(
    unsigned int i) const
{
    return (*points_)[i].cast<scalar_type>();
}

template <class DofType, class PointType, class BasisFunctionType>
inline void
meshless_model_t<DofType, PointType, BasisFunctionType>::in_support_query_t::computeHull(
    unsigned int b,
    unsigned int n,
    Discregrid::BoundingSphere& hull) const
{
    auto vertices_of_sphere = std::vector<Eigen::Vector3d>(n);
    for (unsigned int i = b; i < n + b; ++i)
        vertices_of_sphere[i - b] = (*points_)[m_lst[i]].cast<scalar_type>();

    Discregrid::BoundingSphere const s(vertices_of_sphere);

    hull.x() = s.x();
    hull.r() = s.r() + tolerance_;
}

} // namespace math
} // namespace sbs

#endif // SBS_MATH_MESHLESS_MODEL_H