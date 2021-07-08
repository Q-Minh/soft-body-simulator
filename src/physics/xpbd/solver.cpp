#include "physics/xpbd/solver.h"

#include "physics/collision/intersections.h"
#include "physics/xpbd/collision_constraint.h"
#include "physics/xpbd/distance_constraint.h"
#include "physics/xpbd/green_constraint.h"
#include "physics/xpbd/mesh.h"

#include <numeric>

namespace sbs {
namespace physics {
namespace xpbd {

solver_t::solver_t(double timestep, std::uint32_t substeps, std::uint32_t iterations)
    : physics_bodies_{},
      environment_bodies_{},
      constraints_{},
      collision_constraints_{},
      tetrahedron_to_constraint_map_{},
      dt_{timestep},
      substeps_{substeps},
      iteration_count_{iterations}
{
}

solver_t::solver_t(std::vector<std::shared_ptr<common::renderable_node_t>> const& bodies)
    : physics_bodies_{},
      environment_bodies_{},
      constraints_{},
      collision_constraints_{},
      tetrahedron_to_constraint_map_{},
      dt_{0.0167},
      substeps_{30u},
      iteration_count_{30u}
{
    setup(bodies);
}

solver_t::solver_t(
    double timestep,
    std::uint32_t substeps,
    std::uint32_t iterations,
    std::vector<std::shared_ptr<common::renderable_node_t>> const& bodies)
    : solver_t{timestep, substeps, iterations}
{
    setup(bodies);
}

void solver_t::setup(std::vector<std::shared_ptr<common::renderable_node_t>> const& bodies)
{
    reset();

    for (auto const& body : bodies)
    {
        std::shared_ptr<xpbd::tetrahedral_mesh_t> const physics_body_ptr =
            std::dynamic_pointer_cast<xpbd::tetrahedral_mesh_t>(body);
        if (static_cast<bool>(physics_body_ptr))
        {
            physics_bodies_.push_back(physics_body_ptr);
            continue;
        }

        std::shared_ptr<common::shared_vertex_surface_mesh_i> const environment_body_ptr =
            std::dynamic_pointer_cast<common::shared_vertex_surface_mesh_i>(body);
        if (static_cast<bool>(environment_body_ptr))
        {
            environment_bodies_.push_back(environment_body_ptr);
            continue;
        }
    }

    for (std::shared_ptr<xpbd::tetrahedral_mesh_t> physics_body : physics_bodies_)
    {
        auto const constraint_type = physics_body->simulation_parameters().constraint_type;

        if (constraint_type == constraint_type_t::green)
        {
            create_green_constraints_for_body(physics_body.get());
        }
        if (constraint_type == constraint_type_t::distance)
        {
            create_distance_constraints_for_body(physics_body.get());
        }
    }

    x0_.resize(physics_bodies_.size());
    lagrange_multipliers_.resize(constraints_.size());
}

void solver_t::reset()
{
    physics_bodies_.clear();
    environment_bodies_.clear();
    constraints_.clear();
    collision_constraints_.clear();
    tetrahedron_to_constraint_map_.clear();
    x0_.clear();
    lagrange_multipliers_.clear();
}

void solver_t::step()
{
    // TODO: cut here
    // ...

    bool can_perform_substepping = (iteration_count_ >= substeps_);
    std::uint32_t const num_iterations =
        can_perform_substepping ?
            static_cast<std::uint32_t>(
                static_cast<double>(iteration_count_) / static_cast<double>(substeps_)) :
            iteration_count_;
    std::uint32_t const num_substeps = can_perform_substepping ? substeps_ : 1u;
    double const dt = can_perform_substepping ? dt_ / static_cast<double>(substeps_) : dt_;

    // Hold previous positions x0.
    // Using a vector of vector will result in two levels of indirection before accessing any
    // position. The position vectors will be allocated on the heap at potentially different far
    // locations. We can use a smarter allocation system where these vectors will be stored
    // contiguously.
    std::transform(
        physics_bodies_.begin(),
        physics_bodies_.end(),
        x0_.begin(),
        [](std::shared_ptr<tetrahedral_mesh_t> const& body) {
            std::vector<Eigen::Vector3d> p0{};
            p0.reserve(body->vertices().size());
            std::transform(
                body->vertices().begin(),
                body->vertices().end(),
                std::back_inserter(p0),
                [](vertex_t const& v) { return v.position(); });
            return p0;
        });

    // explicit integration step
    Eigen::Vector3d const gravity{0., -9.81, 0.};
    for (auto const& body : physics_bodies_)
    {
        std::for_each(
            body->vertices().begin(),
            body->vertices().end(),
            [dt, gravity](vertex_t& vertex) {
                vertex.force() += gravity;
                Eigen::Vector3d const acceleration = vertex.force().array() / vertex.mass();
                Eigen::Vector3d const velocity     = vertex.velocity() + dt * acceleration;
                vertex.position() += dt * velocity;
            });
    }

    // TODO: detect collisions
    // ...

    // constraint projection
    std::size_t const J  = constraints_.size();
    std::size_t const Mc = collision_constraints_.size();
    for (std::uint32_t s = 0u; s < num_substeps; ++s)
    {
        std::fill(lagrange_multipliers_.begin(), lagrange_multipliers_.end(), 0.);
        for (std::uint32_t n = 0u; n < num_iterations; ++n)
        {
            for (std::size_t j = 0u; j < J; ++j)
            {
                std::unique_ptr<constraint_t> const& constraint = constraints_[j];
                constraint->project(physics_bodies_, lagrange_multipliers_[j], dt);
            }
            for (std::size_t c = 0u; c < Mc; ++c)
            {
                double infinite_stiffness_lagrange_multiplier   = 0.;
                std::unique_ptr<constraint_t> const& constraint = collision_constraints_[c];
                constraint->project(physics_bodies_, infinite_stiffness_lagrange_multiplier, dt);
            }
        }
    }

    // update positions and velocities
    for (std::size_t b = 0u; b < physics_bodies_.size(); ++b)
    {
        auto const& body               = physics_bodies_[b];
        std::size_t const vertex_count = body->vertices().size();
        for (std::size_t vi = 0u; vi < vertex_count; ++vi)
        {
            vertex_t& vertex         = body->vertices().at(vi);
            Eigen::Vector3d const x0 = x0_[b][vi];
            Eigen::Vector3d const& x = vertex.position();

            vertex.velocity() = (x - x0) / dt;
            vertex.position() = x;
        }
    }

    // friction or other non-conservative forces here

    // reset forces
    for (auto const& body : physics_bodies_)
    {
        for (vertex_t& vertex : body->vertices())
        {
            vertex.force().setZero();
        }
    }
}

double const& solver_t::timestep() const
{
    return dt_;
}

double& solver_t::timestep()
{
    return dt_;
}

std::uint32_t const& solver_t::iterations() const
{
    return iteration_count_;
}

std::uint32_t& solver_t::iterations()
{
    return iteration_count_;
}

std::uint32_t const& solver_t::substeps() const
{
    return substeps_;
}

std::uint32_t& solver_t::substeps()
{
    return substeps_;
}

std::vector<std::shared_ptr<xpbd::tetrahedral_mesh_t>> const& solver_t::simulated_bodies() const
{
    return physics_bodies_;
}

void solver_t::create_green_constraints_for_body(xpbd::tetrahedral_mesh_t* body)
{
    simulation_parameters_t const& params = body->simulation_parameters();
    for (tetrahedron_t const& tetrahedron : body->tetrahedra())
    {
        auto constraint = std::make_unique<green_constraint_t>(
            params.alpha,
            std::make_pair(body, tetrahedron.v1()),
            std::make_pair(body, tetrahedron.v2()),
            std::make_pair(body, tetrahedron.v3()),
            std::make_pair(body, tetrahedron.v4()),
            params.young_modulus,
            params.poisson_ratio);

        constraints_.push_back(std::move(constraint));
    }
}

void solver_t::create_distance_constraints_for_body(xpbd::tetrahedral_mesh_t* body)
{
    simulation_parameters_t const& params = body->simulation_parameters();
    for (edge_t const& edge : body->edges())
    {
        auto constraint = std::make_unique<distance_constraint_t>(
            1. / params.hooke_coefficient,
            std::make_pair(body, edge.v1()),
            std::make_pair(body, edge.v2()));

        constraints_.push_back(std::move(constraint));
    }
}

// void solver_t::step(double timestep, std::uint32_t iterations, std::uint32_t substeps)
//{
//    auto const num_iterations = iterations / substeps;
//    double dt                 = timestep / static_cast<double>(substeps);
//    auto const J              = constraints_.size();
//    std::vector<double> lagrange_multipliers(J, 0.);
//    std::vector<double> collision_lagrange_multipliers{};
//
//    std::vector<Eigen::Matrix3Xd> P{};
//    P.reserve(bodies_->size());
//
//    std::vector<Eigen::VectorXd> M{};
//    M.reserve(bodies_->size());
//
//    for (auto s = 0u; s < substeps; ++s)
//    {
//        P.clear();
//        M.clear();
//
//        /**
//         * Explicit euler step
//         */
//        for (auto const& body : *bodies_)
//        {
//            auto& v = body->physical_model.velocities();
//            auto& x = body->physical_model.positions();
//
//            Eigen::VectorXd const& m     = body->physical_model.masses();
//            Eigen::Matrix3Xd const& fext = body->physical_model.forces();
//            Eigen::Matrix3Xd const a     = fext.array().rowwise() / m.transpose().array();
//
//            auto vexplicit           = v + dt * a;
//            Eigen::Matrix3Xd const p = x + dt * vexplicit;
//            P.push_back(p);
//            M.push_back(m);
//        }
//
//        // generate collision constraints here ...
//        // handle_collisions(P);
//        // std::size_t const Mc = collision_constraints_.size();
//        // collision_lagrange_multipliers.resize(Mc);
//        // std::fill(
//        //    collision_lagrange_multipliers.begin(),
//        //    collision_lagrange_multipliers.end(),
//        //    0.0);
//
//        // sequential gauss seidel type solve
//        std::fill(lagrange_multipliers.begin(), lagrange_multipliers.end(), 0.0);
//        for (auto n = 0u; n < num_iterations; ++n)
//        {
//            for (auto j = 0u; j < J; ++j)
//            {
//                auto const& constraint = constraints_[j];
//                constraint->project(P, M, lagrange_multipliers[j], dt);
//            }
//            // for (auto c = 0u; c < Mc; ++c)
//            //{
//            //    auto const& constraint = collision_constraints_[c];
//            //    constraint->project(P, M, collision_lagrange_multipliers[c], dt);
//            //}
//        }
//
//        // collision_constraints_.clear();
//
//        // update solution
//        for (std::size_t b = 0u; b < bodies_->size(); ++b)
//        {
//            auto const& body = (*bodies_)[b];
//            auto const& p    = P[b];
//            auto const& x    = body->physical_model.positions();
//
//            body->physical_model.velocities() = (p - x) / dt;
//            body->physical_model.positions()  = p;
//        }
//
//        // friction or other non-conservative forces here ...
//    }
//}

// void solver_t::handle_collisions(std::vector<Eigen::Matrix3Xd> const& P)
//{
//    /**
//     * Brute-force collision detection... To be changed
//     * in the future
//     */
//    std::vector<collision::line_segment_t> line_segments{};
//    for (std::size_t bi = 0u; bi < bodies_->size(); ++bi)
//    {
//        auto const& body1 = (*bodies_)[bi];
//        /**
//         * Instead of looking at all positions, should we only
//         * detect collisions for the boundary vertices?
//         */
//        auto const& p = P[bi];
//        auto const& x = body1->physical_model.positions();
//
//        /**
//         * Find line segments x(n) -> p(n+1)
//         */
//        line_segments.clear();
//        line_segments.reserve(p.cols());
//        std::size_t const num_vertices = static_cast<std::size_t>(x.cols());
//        for (std::size_t vi = 0u; vi < num_vertices; ++vi)
//        {
//            line_segments.push_back(collision::line_segment_t{x.col(vi), p.col(vi)});
//        }
//
//        /**
//         * Detect and handle collisions against other physically simulated bodies
//         */
//        for (std::size_t bj = 0u; bj < bodies_->size(); ++bj)
//        {
//            if (bj == bi)
//                continue;
//
//            auto const& body2 = (*bodies_)[bj];
//
//            /**
//             * Handle collisions between line segments and boundary faces of other bodies
//             */
//            auto const& boundary        = body2->render_model.triangles();
//            auto const& index_map       = body2->render_model.index_map();
//            std::size_t const num_faces = static_cast<std::size_t>(boundary.cols());
//            for (std::size_t f = 0u; f < num_faces; ++f)
//            {
//                auto const v1 = index_map[boundary(0u, f)];
//                auto const v2 = index_map[boundary(1u, f)];
//                auto const v3 = index_map[boundary(2u, f)];
//
//                collision::triangle_t const triangle{
//                    body2->physical_model.positions().col(v1),
//                    body2->physical_model.positions().col(v2),
//                    body2->physical_model.positions().col(v3)};
//
//                collision::normal_t const normal = triangle.normal();
//
//                for (std::size_t vi = 0u; vi < num_vertices; ++vi)
//                {
//                    collision::line_segment_t const& line_segment = line_segments[vi];
//
//                    auto const intersection = collision::intersect(line_segment, triangle);
//                    if (!intersection.has_value())
//                        continue;
//
//                    // double const alpha        = per_body_simulation_parameters_[bi].alpha;
//                    auto collision_constraint = std::make_unique<collision_constraint_t>(
//                        // alpha,
//                        std::make_pair(vi, bi),
//                        intersection.value(),
//                        normal);
//
//                    collision_constraints_.push_back(std::move(collision_constraint));
//                }
//            }
//        }
//    }
//}

} // namespace xpbd
} // namespace physics
} // namespace sbs