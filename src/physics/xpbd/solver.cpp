#include "physics/xpbd/solver.h"

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
        for (auto& body : bodies_)
        {
            auto& v = body->mesh.velocities();
            auto& x = body->mesh.positions();

            Eigen::VectorXd const& m     = body->mesh.masses();
            Eigen::Matrix3Xd const& fext = body->mesh.forces();
            Eigen::Matrix3Xd const a = fext.array().rowwise() / m.transpose().array();

            auto vexplicit           = v + dt * a;
            Eigen::Matrix3Xd const p = x + dt * vexplicit;
            P.push_back(std::move(p));
            M.push_back(std::move(m));
        }

        // generate collision constraints here ...
        handle_collisions();

        // sequential gauss seidel type solve
        std::fill(lagrange_multipliers.begin(), lagrange_multipliers.end(), 0.0);
        for (auto n = 0u; n < num_iterations; ++n)
        {
            for (auto j = 0u; j < J; ++j)
            {
                auto const& constraint = constraints_[j];
                constraint->project(P, M, lagrange_multipliers[j], dt);
            }
        }

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
        }

        // friction or other non-conservative forces here ...
    }
}

void solver_t::handle_collisions() {}

} // namespace xpbd
} // namespace physics
} // namespace sbs