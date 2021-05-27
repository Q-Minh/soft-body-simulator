#include "physics/xpbd/solver.h"

#include "physics/collision/intersections.h"
#include "physics/xpbd/collision_constraint.h"
#include "physics/xpbd/distance_constraint.h"
#include "physics/xpbd/green_constraint.h"

#include <numeric>

namespace sbs {
namespace physics {
namespace xpbd {

void solver_t::setup(
    std::vector<std::shared_ptr<common::node_t>> const& bodies,
    std::vector<simulation_parameters_t> const& per_body_simulation_parameters)
{
    bodies_                         = bodies;
    per_body_simulation_parameters_ = per_body_simulation_parameters;

    bool const bodies_and_constraint_types_vectors_are_of_same_size =
        bodies.size() == per_body_simulation_parameters.size();
    assert(bodies_and_constraint_types_vectors_are_of_same_size);

    std::size_t const num_elements = std::accumulate(
        bodies.begin(),
        bodies.end(),
        std::size_t{0u},
        [](std::size_t sum, std::shared_ptr<common::node_t> const& body) {
            std::size_t const num_elements_of_body =
                static_cast<std::size_t>(body->mesh.elements().cols());
            return sum + num_elements_of_body;
        });

    /**
     * Pre-allocate an estimated number of constraints to reduce memory overhead of
     * vector reallocations on resize.
     */
    constraints_.reserve(num_elements);

    /**
     * Initialize constraints
     */
    for (std::size_t b = 0u; b < bodies.size(); ++b)
    {
        std::shared_ptr<common::node_t> const& body = bodies_[b];

        if (body->is_fixed)
            continue;

        constraint_type_t const constraint_type =
            per_body_simulation_parameters_[b].constraint_type;

        std::size_t const num_elements_of_body =
            static_cast<std::size_t>(body->mesh.elements().cols());

        double const alpha = per_body_simulation_parameters_[b].alpha;

        if (constraint_type == constraint_type_t::distance)
        {
            auto const edges = common::edges(body->mesh);
            std::transform(
                edges.begin(),
                edges.end(),
                std::back_inserter(constraints_),
                [alpha, b, body](std::pair<std::uint32_t, std::uint32_t> const& edge) {
                    return std::make_unique<distance_constraint_t>(
                        alpha,
                        body->mesh.positions(),
                        std::make_pair(edge.first, b),
                        std::make_pair(edge.second, b));
                });
        }
        if (constraint_type == constraint_type_t::green)
        {
            for (std::size_t e = 0u; e < num_elements_of_body; ++e)
            {
                auto const v1 = body->mesh.elements()(0u, e);
                auto const v2 = body->mesh.elements()(1u, e);
                auto const v3 = body->mesh.elements()(2u, e);
                auto const v4 = body->mesh.elements()(3u, e);

                double const young_modulus = per_body_simulation_parameters_[b].young_modulus;
                double const poisson_ratio = per_body_simulation_parameters_[b].poisson_ratio;

                auto constraint = std::make_unique<green_constraint_t>(
                    alpha,
                    body->mesh.positions(),
                    std::make_pair(v1, b),
                    std::make_pair(v2, b),
                    std::make_pair(v3, b),
                    std::make_pair(v4, b),
                    young_modulus,
                    poisson_ratio);

                constraints_.push_back(std::move(constraint));
            }
        }
    }
}

void solver_t::step(double timestep, std::uint32_t iterations, std::uint32_t substeps)
{
    auto const num_iterations = iterations / substeps;
    double dt                 = timestep / static_cast<double>(substeps);
    auto const J              = constraints_.size();
    std::vector<double> lagrange_multipliers(J, 0.);
    std::vector<double> collision_lagrange_multipliers{};

    std::vector<Eigen::Matrix3Xd> P{};
    P.reserve(bodies_.size());

    std::vector<Eigen::VectorXd> M{};
    M.reserve(bodies_.size());

    for (auto s = 0u; s < substeps; ++s)
    {
        P.clear();
        M.clear();

        /**
         * Explicit euler step
         */
        for (auto const& body : bodies_)
        {
            auto& v = body->mesh.velocities();
            auto& x = body->mesh.positions();

            Eigen::VectorXd const& m     = body->mesh.masses();
            Eigen::Matrix3Xd const& fext = body->mesh.forces();
            Eigen::Matrix3Xd const a     = fext.array().rowwise() / m.transpose().array();

            auto vexplicit           = v + dt * a;
            Eigen::Matrix3Xd const p = x + dt * vexplicit;
            P.push_back(std::move(p));
            M.push_back(std::move(m));
        }

        // generate collision constraints here ...
        handle_collisions(P);
        std::size_t const Mc = collision_constraints_.size();
        collision_lagrange_multipliers.resize(Mc);
        std::fill(
            collision_lagrange_multipliers.begin(),
            collision_lagrange_multipliers.end(),
            0.0);

        // sequential gauss seidel type solve
        std::fill(lagrange_multipliers.begin(), lagrange_multipliers.end(), 0.0);
        for (auto n = 0u; n < num_iterations; ++n)
        {
            for (auto j = 0u; j < J; ++j)
            {
                auto const& constraint = constraints_[j];
                constraint->project(P, M, lagrange_multipliers[j], dt);
            }
            for (auto c = 0u; c < Mc; ++c)
            {
                auto const& constraint = collision_constraints_[c];
                constraint->project(P, M, collision_lagrange_multipliers[c], dt);
            }
        }

        collision_constraints_.clear();

        // update solution
        for (std::size_t b = 0u; b < bodies_.size(); ++b)
        {
            auto const& body = bodies_[b];
            auto const& p    = P[b];
            auto const& x    = body->mesh.positions();

            if (body->is_fixed)
                continue;

            body->mesh.velocities() = (p - x) / dt;
            body->mesh.positions()  = p;
            body->render_state      = sbs::common::node_t::render_state_t::dirty;
        }

        // friction or other non-conservative forces here ...
    }
}

void solver_t::handle_collisions(std::vector<Eigen::Matrix3Xd> const& P)
{
    /**
     * Brute-force collision detection... To be changed
     * in the future
     */
    std::vector<collision::line_segment_t> line_segments{};
    for (std::size_t bi = 0u; bi < bodies_.size(); ++bi)
    {
        auto const& body1 = bodies_[bi];
        /**
         * Instead of looking at all positions, should we only
         * detect collisions for the boundary vertices?
         */
        auto const& p = P[bi];
        auto const& x = body1->mesh.positions();

        /**
         * Find line segments x(n) -> p(n+1)
         */
        line_segments.clear();
        line_segments.reserve(p.cols());
        std::size_t const num_vertices = static_cast<std::size_t>(x.cols());
        for (std::size_t vi = 0u; vi < num_vertices; ++vi)
        {
            line_segments.push_back(collision::line_segment_t{x.col(vi), p.col(vi)});
        }

        for (std::size_t bj = 0u; bj < bodies_.size(); ++bj)
        {
            if (bj == bi)
                continue;

            auto const& body2 = bodies_[bj];

            /**
             * Handle collisions between line segments and boundary faces of other bodies
             */
            auto const& boundary        = body2->mesh.faces();
            std::size_t const num_faces = static_cast<std::size_t>(boundary.cols());
            for (std::size_t f = 0u; f < num_faces; ++f)
            {
                auto const v1 = boundary(0u, f);
                auto const v2 = boundary(1u, f);
                auto const v3 = boundary(2u, f);

                collision::triangle_t const triangle{
                    body2->mesh.positions().col(v1),
                    body2->mesh.positions().col(v2),
                    body2->mesh.positions().col(v3)};

                collision::normal_t const normal = triangle.normal();

                for (std::size_t vi = 0u; vi < num_vertices; ++vi)
                {
                    collision::line_segment_t const& line_segment = line_segments[vi];

                    auto const intersection = collision::intersect(line_segment, triangle);
                    if (!intersection.has_value())
                        continue;

                    // double const alpha        = per_body_simulation_parameters_[bi].alpha;
                    auto collision_constraint = std::make_unique<collision_constraint_t>(
                        // alpha,
                        std::make_pair(vi, bi),
                        intersection.value(),
                        normal);

                    collision_constraints_.push_back(std::move(collision_constraint));
                }
            }
        }
    }
}

} // namespace xpbd
} // namespace physics
} // namespace sbs