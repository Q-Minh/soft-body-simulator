#include "physics/cutting/cut_tetrahedron.h"

#include "common/primitive.h"

#include <iostream>

namespace sbs {
namespace physics {
namespace cutting {

static double lerp_coefficient(
    Eigen::Vector3d const& A,
    Eigen::Vector3d const& B,
    Eigen::Vector3d const& P,
    double const eps = 1e-8)
{
    double const dx = B.x() - A.x();
    if (std::abs(dx) > eps)
        return (P.x() - A.x()) / dx;

    double const dy = B.y() - A.y();
    if (std::abs(dy) > eps)
        return (P.y() - A.y()) / dy;

    double const dz = B.z() - A.z();
    if (std::abs(dz) > eps)
        return (P.z() - A.z()) / dz;

    return 0.;
}

tetrahedral_mesh_cutter_t::tetrahedral_mesh_cutter_t(common::shared_vertex_mesh_t& mesh)
    : mesh_(mesh),
      cut_edges_(),
      cut_faces_(),
      cut_tetrahedra_(),
      previous_vertex_count_(static_cast<std::uint32_t>(mesh.positions().cols()))
{
}

void tetrahedral_mesh_cutter_t::detect_cuts(
    common::shared_vertex_surface_mesh_t const& cutting_surface)
{
    std::uint32_t const num_tetrahedra = static_cast<std::uint32_t>(mesh_.elements().cols());
    for (std::uint32_t e = 0u; e < num_tetrahedra; ++e)
    {
        detect_cuts(e, cutting_surface);
    }
}

void tetrahedral_mesh_cutter_t::detect_cuts(
    std::uint32_t tetrahedron,
    common::shared_vertex_surface_mesh_t const& cutting_surface)
{
    auto const v1 = mesh_.elements()(0u, tetrahedron);
    auto const v2 = mesh_.elements()(1u, tetrahedron);
    auto const v3 = mesh_.elements()(2u, tetrahedron);
    auto const v4 = mesh_.elements()(3u, tetrahedron);

    auto const& p1 = mesh_.positions().col(v1);
    auto const& p2 = mesh_.positions().col(v2);
    auto const& p3 = mesh_.positions().col(v3);
    auto const& p4 = mesh_.positions().col(v4);

    collision::line_segment_t const e1{p1, p2};
    collision::line_segment_t const e2{p2, p3};
    collision::line_segment_t const e3{p3, p1};
    collision::line_segment_t const e4{p1, p4};
    collision::line_segment_t const e5{p2, p4};
    collision::line_segment_t const e6{p3, p4};

    std::byte edge_intersection_mask{0b00000000};

    /**
     * Find tetrahedron's cut edges
     */
    std::size_t const num_cutting_triangles = cutting_surface.triangle_count();

    std::size_t const num_vertices = previous_vertex_count_;

    auto const get_num_new_vertices = [this]() {
        return 2u * cut_edges_.size() + cut_faces_.size();
    };

    auto const is_over_triangle = [](Eigen::Vector3d const& p, collision::triangle_t const& t) {
        auto const n = t.normal();
        auto const v = p - t.a; // Take p from any point on the plane spanned by t
        return v.dot(n) > 0.;
    };

    auto const create_cut_edge_facet = [num_vertices, get_num_new_vertices, is_over_triangle, this](
                                           std::uint32_t ev1,
                                           std::uint32_t ev2,
                                           Eigen::Vector3d const& intersection,
                                           collision::triangle_t const& cut) {
        edge_facet_key_type const edge_key{ev1, ev2};
        bool const is_first_detection = (cut_edges_.find(edge_key) == cut_edges_.end());
        if (is_first_detection)
        {
            std::uint32_t const new_vertex_index =
                static_cast<std::uint32_t>(num_vertices + get_num_new_vertices());
            cut_edges_[edge_key].vi = new_vertex_index;

            auto const& ep1                 = mesh_.positions().col(ev1);
            bool const is_ev1_over_triangle = is_over_triangle(ep1, cut);

            cut_edges_[edge_key].v1           = is_ev1_over_triangle ? ev2 : ev1;
            cut_edges_[edge_key].v2           = is_ev1_over_triangle ? ev1 : ev2;
            cut_edges_[edge_key].intersection = intersection;
        }
    };

    for (std::size_t f = 0u; f < num_cutting_triangles; ++f)
    {
        auto const fv1 = cutting_surface.triangles()(0u, f);
        auto const fv2 = cutting_surface.triangles()(1u, f);
        auto const fv3 = cutting_surface.triangles()(2u, f);

        collision::triangle_t const cutting_triangle{
            cutting_surface.vertices().col(fv1),
            cutting_surface.vertices().col(fv2),
            cutting_surface.vertices().col(fv3)};

        auto const e1_intersection = collision::intersect_twoway(e1, cutting_triangle);
        auto const e2_intersection = collision::intersect_twoway(e2, cutting_triangle);
        auto const e3_intersection = collision::intersect_twoway(e3, cutting_triangle);
        auto const e4_intersection = collision::intersect_twoway(e4, cutting_triangle);
        auto const e5_intersection = collision::intersect_twoway(e5, cutting_triangle);
        auto const e6_intersection = collision::intersect_twoway(e6, cutting_triangle);

        if (e1_intersection.has_value())
        {
            edge_intersection_mask |= std::byte{0b00000001};
            create_cut_edge_facet(v1, v2, *e1_intersection, cutting_triangle);
        }
        if (e2_intersection.has_value())
        {
            edge_intersection_mask |= std::byte{0b00000010};
            create_cut_edge_facet(v2, v3, *e2_intersection, cutting_triangle);
        }
        if (e3_intersection.has_value())
        {
            edge_intersection_mask |= std::byte{0b00000100};
            create_cut_edge_facet(v3, v1, *e3_intersection, cutting_triangle);
        }
        if (e4_intersection.has_value())
        {
            edge_intersection_mask |= std::byte{0b00001000};
            create_cut_edge_facet(v1, v4, *e4_intersection, cutting_triangle);
        }
        if (e5_intersection.has_value())
        {
            edge_intersection_mask |= std::byte{0b00010000};
            create_cut_edge_facet(v2, v4, *e5_intersection, cutting_triangle);
        }
        if (e6_intersection.has_value())
        {
            edge_intersection_mask |= std::byte{0b00100000};
            create_cut_edge_facet(v3, v4, *e6_intersection, cutting_triangle);
        }
    }

    collision::triangle_t const f1{p1, p2, p4};
    collision::triangle_t const f2{p2, p3, p4};
    collision::triangle_t const f3{p3, p1, p4};
    collision::triangle_t const f4{p1, p3, p2};

    auto const boundary_edges = cutting_surface.boundary_edges();
    for (auto const& [fv1, fv2] : boundary_edges)
    {
        auto const create_cut_face_facet = [this, num_vertices, get_num_new_vertices](
                                               std::uint32_t fv1,
                                               std::uint32_t fv2,
                                               std::uint32_t fv3,
                                               Eigen::Vector3d const& intersection) {
            triangle_facet_key_type const face_key{fv1, fv2, fv3};
            bool const is_first_detection = (cut_faces_.find(face_key) == cut_faces_.end());
            if (is_first_detection)
            {
                std::uint32_t const new_vertex_index =
                    static_cast<std::uint32_t>(num_vertices + get_num_new_vertices());
                cut_faces_[face_key].vi           = new_vertex_index;
                cut_faces_[face_key].intersection = intersection;
            }
        };

        collision::line_segment_t const boundary_edge{
            cutting_surface.vertices().col(fv1),
            cutting_surface.vertices().col(fv2)};

        auto const f1_intersection = collision::intersect_twoway(boundary_edge, f1);
        auto const f2_intersection = collision::intersect_twoway(boundary_edge, f2);
        auto const f3_intersection = collision::intersect_twoway(boundary_edge, f3);
        auto const f4_intersection = collision::intersect_twoway(boundary_edge, f4);

        if (f1_intersection.has_value())
        {
            create_cut_face_facet(v1, v2, v4, *f1_intersection);
        }
        if (f2_intersection.has_value())
        {
            create_cut_face_facet(v2, v3, v4, *f2_intersection);
        }
        if (f3_intersection.has_value())
        {
            create_cut_face_facet(v3, v1, v4, *f3_intersection);
        }
        if (f4_intersection.has_value())
        {
            create_cut_face_facet(v1, v3, v2, *f4_intersection);
        }
    }

    bool const should_cut_tet = (edge_intersection_mask != std::byte{0b00000000});
    if (!should_cut_tet)
        return;

    tetrahedron_key_type const tetrahedron_key = tetrahedron;
    // cut_tetrahedra_[tetrahedron_key].edge_intersection_mask |= edge_intersection_mask;
    // cut_tetrahedra_[tetrahedron_key].v1 = mesh_.elements()(0u, tetrahedron);
    // cut_tetrahedra_[tetrahedron_key].v2 = mesh_.elements()(1u, tetrahedron);
    // cut_tetrahedra_[tetrahedron_key].v3 = mesh_.elements()(2u, tetrahedron);
    // cut_tetrahedra_[tetrahedron_key].v4 = mesh_.elements()(3u, tetrahedron);
    cut_tetrahedra_[tetrahedron_key] |= edge_intersection_mask;
}

void tetrahedral_mesh_cutter_t::create_new_vertices()
{
    std::uint32_t const num_vertices =
        previous_vertex_count_ +
        static_cast<std::uint32_t>(2u * cut_edges_.size() + cut_faces_.size());

    mesh_.positions().conservativeResize(3u, num_vertices);
    mesh_.forces().conservativeResize(3u, num_vertices);
    mesh_.velocities().conservativeResize(3u, num_vertices);
    mesh_.masses().conservativeResize(num_vertices);

    for (auto& [key, facet] : cut_edges_)
    {
        if (facet.has_created_new_vertex)
            continue;

        auto const& v1 = key[0];
        auto const& v2 = key[1];

        auto const& p1 = mesh_.positions().col(v1);
        auto const& p2 = mesh_.positions().col(v2);

        auto const t = lerp_coefficient(p1, p2, facet.intersection);

        int const vi              = static_cast<int>(facet.vi);
        mesh_.positions().col(vi) = facet.intersection;
        mesh_.forces().col(vi)    = (1 - t) * mesh_.forces().col(v1) + t * mesh_.forces().col(v2);
        mesh_.velocities().col(vi) =
            (1 - t) * mesh_.velocities().col(v1) + t * mesh_.velocities().col(v2);
        mesh_.masses()(vi) = std::min(mesh_.masses()(v1), mesh_.masses()(v2));

        int const vip               = vi + 1;
        mesh_.positions().col(vip)  = mesh_.positions().col(vi);
        mesh_.forces().col(vip)     = mesh_.forces().col(vi);
        mesh_.velocities().col(vip) = mesh_.velocities().col(vi);
        mesh_.masses()(vip)         = mesh_.masses()(vi);

        facet.has_created_new_vertex = true;
    }

    for (auto& [key, facet] : cut_faces_)
    {
        if (facet.has_created_new_vertex)
            continue;

        auto const& v1 = key[0];
        auto const& v2 = key[1];
        auto const& v3 = key[2];

        auto const& A = mesh_.positions().col(v1);
        auto const& B = mesh_.positions().col(v2);
        auto const& C = mesh_.positions().col(v3);

        auto const [u, v, w] = common::barycentric_coordinates(A, B, C, facet.intersection);

        int const vi              = static_cast<int>(facet.vi);
        mesh_.positions().col(vi) = facet.intersection;
        mesh_.forces().col(vi) =
            u * mesh_.forces().col(v1) + v * mesh_.forces().col(v2) + w * mesh_.forces().col(v3);
        mesh_.velocities().col(vi) = u * mesh_.velocities().col(v1) +
                                     v * mesh_.velocities().col(v2) +
                                     w * mesh_.velocities().col(v3);
        mesh_.masses()(vi) =
            u * mesh_.masses()(v1) + v * mesh_.masses()(v2) + w * mesh_.masses()(v3);

        facet.has_created_new_vertex = true;
    }
}

bool tetrahedral_mesh_cutter_t::update_topology(
    common::shared_vertex_surface_mesh_t const& cutting_surface)
{
    bool has_mesh_been_modified{false};

    struct cutting_list_element_type
    {
        std::uint32_t tetrahedron;
        std::uint32_t offset;
        subdivided_element_type subdivided_elements;
    };
    std::list<cutting_list_element_type> cutting_list{};
    std::uint32_t updated_tetrahedra_count = static_cast<std::uint32_t>(mesh_.elements().cols());
    for (auto it = cut_tetrahedra_.begin(); it != cut_tetrahedra_.end();)
    {
        auto const& [tetrahedron, edge_intersection_mask] = *it;

        /**
         * If the cutting surface and this tetrahedron still intersect,
         * we delay the cut. When the cutting surface will leave this
         * tetrahedron in the future, this tetrahedron will be cut.
         */
        if (intersects(tetrahedron, cutting_surface))
        {
            ++it;
            continue;
        }

        auto const tet_to_subdivide = subdivide_mesh(edge_intersection_mask, tetrahedron);
        if (tet_to_subdivide.has_value())
        {
            cutting_list_element_type cutting_list_element{};
            cutting_list_element.tetrahedron         = tetrahedron;
            cutting_list_element.subdivided_elements = *tet_to_subdivide;
            cutting_list_element.offset              = updated_tetrahedra_count;

            std::uint32_t const num_subdivided_elements =
                static_cast<std::uint32_t>(cutting_list_element.subdivided_elements.cols());

            cutting_list.push_back(cutting_list_element);
            /**
             * 1 tetrahedron being subdivided into N elements results in
             * (N - 1) elements to add in our mesh, because we remove 1 tet, and
             * we add back N tets.
             */
            updated_tetrahedra_count += num_subdivided_elements - 1u;
            has_mesh_been_modified = true;
        }

        it = cut_tetrahedra_.erase(it);
    }

    mesh_.elements().conservativeResize(4u, updated_tetrahedra_count);
    for (auto const& cutting_list_element : cutting_list)
    {
        std::uint32_t const tetrahedron     = cutting_list_element.tetrahedron;
        tetrahedra_type const& tetrahedra   = cutting_list_element.subdivided_elements;
        std::uint32_t const offset          = cutting_list_element.offset;
        int const num_additional_tetrahedra = tetrahedra.cols() - 1;

        mesh_.elements().col(tetrahedron) = tetrahedra.col(0u);
        mesh_.elements().block(0, offset, 4, num_additional_tetrahedra) =
            tetrahedra.block(0, 1, 4, num_additional_tetrahedra);
    }

    return has_mesh_been_modified;
}

// std::optional<tetrahedron_mesh_cutter_t::subdivided_element_type>
// tetrahedron_mesh_cutter_t::cut_tetrahedron(
//    std::uint32_t tetrahedron,
//    std::byte const edge_intersection_mask,
//    std::array<Eigen::Vector3d, 6u> const& edge_intersections,
//    std::array<Eigen::Vector3d, 4u> const& face_intersections)
//{
//    auto const tets = subdivide_mesh(
//        edge_intersection_mask,
//        mesh_,
//        tetrahedron,
//        edge_intersections,
//        face_intersections);
//
//    return tets;
//}

// std::optional<tetrahedron_mesh_cutter_t::subdivided_element_type>
// tetrahedron_mesh_cutter_t::cut_tetrahedron(
//    std::uint32_t tetrahedron,
//    common::shared_vertex_surface_mesh_t const& cutting_surface)
//{
//    auto const v1 = mesh_.elements()(0u, tetrahedron);
//    auto const v2 = mesh_.elements()(1u, tetrahedron);
//    auto const v3 = mesh_.elements()(2u, tetrahedron);
//    auto const v4 = mesh_.elements()(3u, tetrahedron);
//
//    auto const& p1 = mesh_.positions().col(v1);
//    auto const& p2 = mesh_.positions().col(v2);
//    auto const& p3 = mesh_.positions().col(v3);
//    auto const& p4 = mesh_.positions().col(v4);
//
//    collision::line_segment_t const e1{p1, p2};
//    collision::line_segment_t const e2{p2, p3};
//    collision::line_segment_t const e3{p3, p1};
//    collision::line_segment_t const e4{p1, p4};
//    collision::line_segment_t const e5{p2, p4};
//    collision::line_segment_t const e6{p3, p4};
//
//    std::byte edge_intersection_mask{0b00000000};
//
//    std::array<Eigen::Vector3d, 6u> edge_intersections{};
//    std::array<Eigen::Vector3d, 4u> face_intersections{};
//
//    /**
//     * Find tetrahedron's cut edges
//     */
//    std::size_t const num_cutting_triangles =
//        static_cast<std::size_t>(cutting_surface.triangles().cols());
//
//    for (std::size_t f = 0u; f < num_cutting_triangles; ++f)
//    {
//        auto const v1 = cutting_surface.triangles()(0u, f);
//        auto const v2 = cutting_surface.triangles()(1u, f);
//        auto const v3 = cutting_surface.triangles()(2u, f);
//
//        collision::triangle_t const cutting_triangle{
//            cutting_surface.vertices().col(v1),
//            cutting_surface.vertices().col(v2),
//            cutting_surface.vertices().col(v3)};
//
//        auto const e1_intersection = collision::intersect_twoway(e1, cutting_triangle);
//        auto const e2_intersection = collision::intersect_twoway(e2, cutting_triangle);
//        auto const e3_intersection = collision::intersect_twoway(e3, cutting_triangle);
//        auto const e4_intersection = collision::intersect_twoway(e4, cutting_triangle);
//        auto const e5_intersection = collision::intersect_twoway(e5, cutting_triangle);
//        auto const e6_intersection = collision::intersect_twoway(e6, cutting_triangle);
//
//        if (e1_intersection.has_value())
//        {
//            edge_intersections[0] = e1_intersection.value();
//            edge_intersection_mask |= std::byte{0b00000001};
//        }
//        if (e2_intersection.has_value())
//        {
//            edge_intersections[1] = e2_intersection.value();
//            edge_intersection_mask |= std::byte{0b00000010};
//        }
//        if (e3_intersection.has_value())
//        {
//            edge_intersections[2] = e3_intersection.value();
//            edge_intersection_mask |= std::byte{0b00000100};
//        }
//        if (e4_intersection.has_value())
//        {
//            edge_intersections[3] = e4_intersection.value();
//            edge_intersection_mask |= std::byte{0b00001000};
//        }
//        if (e5_intersection.has_value())
//        {
//            edge_intersections[4] = e5_intersection.value();
//            edge_intersection_mask |= std::byte{0b00010000};
//        }
//        if (e6_intersection.has_value())
//        {
//            edge_intersections[5] = e6_intersection.value();
//            edge_intersection_mask |= std::byte{0b00100000};
//        }
//    }
//
//    collision::triangle_t const f1{p1, p2, p4};
//    collision::triangle_t const f2{p2, p3, p4};
//    collision::triangle_t const f3{p3, p1, p4};
//    collision::triangle_t const f4{p1, p3, p2};
//
//    auto const boundary_edges = cutting_surface.boundary_edges();
//    for (auto const& [v1, v2] : boundary_edges)
//    {
//        collision::line_segment_t const boundary_edge{
//            cutting_surface.vertices().col(v1),
//            cutting_surface.vertices().col(v2)};
//
//        auto const f1_intersection = collision::intersect_twoway(boundary_edge, f1);
//        auto const f2_intersection = collision::intersect_twoway(boundary_edge, f2);
//        auto const f3_intersection = collision::intersect_twoway(boundary_edge, f3);
//        auto const f4_intersection = collision::intersect_twoway(boundary_edge, f4);
//
//        if (f1_intersection.has_value())
//            face_intersections[0] = f1_intersection.value();
//        if (f2_intersection.has_value())
//            face_intersections[1] = f2_intersection.value();
//        if (f3_intersection.has_value())
//            face_intersections[2] = f3_intersection.value();
//        if (f4_intersection.has_value())
//            face_intersections[3] = f4_intersection.value();
//    }
//
//    return cut_tetrahedron(
//        tetrahedron,
//        edge_intersection_mask,
//        edge_intersections,
//        face_intersections);
//}

bool tetrahedral_mesh_cutter_t ::cut_tetrahedral_mesh(
    common::shared_vertex_surface_mesh_t const& cutting_surface)
{
    /**
     * Brute-force search for intersections between the cutting surface and the mesh.
     * Pruning the search space using efficient spatial acceleration structures
     * could be beneficial for performance.
     */
    detect_cuts(cutting_surface);

    create_new_vertices();

    bool const has_topology_changed = update_topology(cutting_surface);
    // if (!has_topology_changed)
    //    finalize_cut();

    return has_topology_changed;
}

void tetrahedral_mesh_cutter_t::finalize_cut()
{
    cut_edges_.clear();
    cut_faces_.clear();
    cut_tetrahedra_.clear();
    previous_vertex_count_ = static_cast<std::uint32_t>(mesh_.positions().cols());
}

bool tetrahedral_mesh_cutter_t::intersects(
    std::uint32_t tetrahedron,
    common::shared_vertex_surface_mesh_t const& cutting_surface) const
{
    auto const v1 = mesh_.elements()(0u, tetrahedron);
    auto const v2 = mesh_.elements()(1u, tetrahedron);
    auto const v3 = mesh_.elements()(2u, tetrahedron);
    auto const v4 = mesh_.elements()(3u, tetrahedron);

    auto const& p1 = mesh_.positions().col(v1);
    auto const& p2 = mesh_.positions().col(v2);
    auto const& p3 = mesh_.positions().col(v3);
    auto const& p4 = mesh_.positions().col(v4);

    collision::line_segment_t const e1{p1, p2};
    collision::line_segment_t const e2{p2, p3};
    collision::line_segment_t const e3{p3, p1};
    collision::line_segment_t const e4{p1, p4};
    collision::line_segment_t const e5{p2, p4};
    collision::line_segment_t const e6{p3, p4};

    std::size_t const num_cutting_triangles = cutting_surface.triangle_count();

    for (std::size_t f = 0u; f < num_cutting_triangles; ++f)
    {
        auto const fv1 = cutting_surface.triangles()(0u, f);
        auto const fv2 = cutting_surface.triangles()(1u, f);
        auto const fv3 = cutting_surface.triangles()(2u, f);

        collision::triangle_t const cutting_triangle{
            cutting_surface.vertices().col(fv1),
            cutting_surface.vertices().col(fv2),
            cutting_surface.vertices().col(fv3)};

        auto const e1_intersection = collision::intersect_twoway(e1, cutting_triangle);
        if (e1_intersection.has_value())
            return true;

        auto const e2_intersection = collision::intersect_twoway(e2, cutting_triangle);
        if (e2_intersection.has_value())
            return true;

        auto const e3_intersection = collision::intersect_twoway(e3, cutting_triangle);
        if (e3_intersection.has_value())
            return true;

        auto const e4_intersection = collision::intersect_twoway(e4, cutting_triangle);
        if (e4_intersection.has_value())
            return true;

        auto const e5_intersection = collision::intersect_twoway(e5, cutting_triangle);
        if (e5_intersection.has_value())
            return true;

        auto const e6_intersection = collision::intersect_twoway(e6, cutting_triangle);
        if (e6_intersection.has_value())
            return true;
    }

    collision::triangle_t const f1{p1, p2, p4};
    collision::triangle_t const f2{p2, p3, p4};
    collision::triangle_t const f3{p3, p1, p4};
    collision::triangle_t const f4{p1, p3, p2};

    auto const boundary_edges = cutting_surface.boundary_edges();
    for (auto const& [fv1, fv2] : boundary_edges)
    {
        collision::line_segment_t const boundary_edge{
            cutting_surface.vertices().col(fv1),
            cutting_surface.vertices().col(fv2)};

        auto const f1_intersection = collision::intersect_twoway(boundary_edge, f1);
        if (f1_intersection.has_value())
            return true;

        auto const f2_intersection = collision::intersect_twoway(boundary_edge, f2);
        if (f2_intersection.has_value())
            return true;

        auto const f3_intersection = collision::intersect_twoway(boundary_edge, f3);
        if (f3_intersection.has_value())
            return true;

        auto const f4_intersection = collision::intersect_twoway(boundary_edge, f4);
        if (f4_intersection.has_value())
            return true;
    }

    return false;
}

std::optional<tetrahedral_mesh_cutter_t::subdivided_element_type>
tetrahedral_mesh_cutter_t::subdivide_mesh(
    std::byte const& edge_intersection_mask,
    std::uint32_t tetrahedron)
{
    auto const _v1 = mesh_.elements()(0u, tetrahedron);
    auto const _v2 = mesh_.elements()(1u, tetrahedron);
    auto const _v3 = mesh_.elements()(2u, tetrahedron);
    auto const _v4 = mesh_.elements()(3u, tetrahedron);

    edge_facet_key_type const e1_key{_v1, _v2};
    edge_facet_key_type const e2_key{_v2, _v3};
    edge_facet_key_type const e3_key{_v3, _v1};
    edge_facet_key_type const e4_key{_v1, _v4};
    edge_facet_key_type const e5_key{_v2, _v4};
    edge_facet_key_type const e6_key{_v3, _v4};

    triangle_facet_key_type const f1_key{_v1, _v2, _v4};
    triangle_facet_key_type const f2_key{_v2, _v3, _v4};
    triangle_facet_key_type const f3_key{_v3, _v1, _v4};
    triangle_facet_key_type const f4_key{_v1, _v3, _v2};

    std::array<edge_facet_key_type, 6u> const
        cut_edge_keys{e1_key, e2_key, e3_key, e4_key, e5_key, e6_key};
    std::array<triangle_facet_key_type, 4u> const cut_face_keys{f1_key, f2_key, f3_key, f4_key};

    int constexpr v1{0};
    int constexpr v2{1};
    int constexpr v3{2};
    int constexpr v4{3};

    int constexpr e1{0};
    int constexpr e2{1};
    int constexpr e3{2};
    int constexpr e4{3};
    int constexpr e5{4};
    int constexpr e6{5};

    int constexpr f1{0};
    int constexpr f2{1};
    int constexpr f3{2};
    int constexpr f4{3};

    std::byte constexpr case_1_125{0b00010011};
    std::byte constexpr case_1_134{0b00001101};
    std::byte constexpr case_1_236{0b00100110};
    std::byte constexpr case_1_456{0b00111000};

    if (edge_intersection_mask == case_1_125)
    {
        auto const& ce1 = cut_edge_keys[e1];
        auto const& ce2 = cut_edge_keys[e2];
        auto const& ce3 = cut_edge_keys[e5];
        auto const tets =
            subdivide_mesh_for_common_case_1(tetrahedron, {v1, v3, v4, v2}, {ce1, ce2, ce3});
        return tets;
    }
    if (edge_intersection_mask == case_1_134)
    {
        auto const& ce1 = cut_edge_keys[e1];
        auto const& ce2 = cut_edge_keys[e4];
        auto const& ce3 = cut_edge_keys[e3];
        auto const tets =
            subdivide_mesh_for_common_case_1(tetrahedron, {v2, v4, v3, v1}, {ce1, ce2, ce3});
        return tets;
    }
    if (edge_intersection_mask == case_1_236)
    {
        auto const& ce1 = cut_edge_keys[e2];
        auto const& ce2 = cut_edge_keys[e3];
        auto const& ce3 = cut_edge_keys[e6];
        auto const tets =
            subdivide_mesh_for_common_case_1(tetrahedron, {v2, v1, v4, v3}, {ce1, ce2, ce3});
        return tets;
    }
    if (edge_intersection_mask == case_1_456)
    {
        auto const& ce1 = cut_edge_keys[e4];
        auto const& ce2 = cut_edge_keys[e5];
        auto const& ce3 = cut_edge_keys[e6];
        auto const tets =
            subdivide_mesh_for_common_case_1(tetrahedron, {v1, v2, v3, v4}, {ce1, ce2, ce3});
        return tets;
    }

    std::byte constexpr case_2_1246{0b00101011};
    std::byte constexpr case_2_1356{0b00110101};
    std::byte constexpr case_2_2345{0b00011110};

    if (edge_intersection_mask == case_2_1246)
    {
        auto const& ce1 = cut_edge_keys[e1];
        auto const& ce2 = cut_edge_keys[e2];
        auto const& ce3 = cut_edge_keys[e6];
        auto const& ce4 = cut_edge_keys[e4];
        auto const tets =
            subdivide_mesh_for_common_case_2(tetrahedron, {v2, v4, v3, v1}, {ce1, ce2, ce3, ce4});
        return tets;
    }
    if (edge_intersection_mask == case_2_1356)
    {
        auto const& ce1 = cut_edge_keys[e5];
        auto const& ce2 = cut_edge_keys[e1];
        auto const& ce3 = cut_edge_keys[e3];
        auto const& ce4 = cut_edge_keys[e6];
        auto const tets =
            subdivide_mesh_for_common_case_2(tetrahedron, {v2, v3, v1, v4}, {ce1, ce2, ce3, ce4});
        return tets;
    }
    if (edge_intersection_mask == case_2_2345)
    {
        auto const& ce1 = cut_edge_keys[e4];
        auto const& ce2 = cut_edge_keys[e3];
        auto const& ce3 = cut_edge_keys[e2];
        auto const& ce4 = cut_edge_keys[e5];
        auto const tets =
            subdivide_mesh_for_common_case_2(tetrahedron, {v1, v2, v3, v4}, {ce1, ce2, ce3, ce4});
        return tets;
    }

    std::byte constexpr case_3_1{0b00000001};
    std::byte constexpr case_3_2{0b00000010};
    std::byte constexpr case_3_3{0b00000100};
    std::byte constexpr case_3_4{0b00001000};
    std::byte constexpr case_3_5{0b00010000};
    std::byte constexpr case_3_6{0b00100000};

    if (edge_intersection_mask == case_3_1)
    {
        auto const& ce1 = cut_edge_keys[e1];
        auto const& cf1 = cut_face_keys[f4];
        auto const& cf2 = cut_face_keys[f1];
        auto const tets =
            subdivide_mesh_for_common_case_3(tetrahedron, {v1, v3, v4, v2}, {ce1}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_3_2)
    {
        auto const& ce1 = cut_edge_keys[e2];
        auto const& cf1 = cut_face_keys[f2];
        auto const& cf2 = cut_face_keys[f4];
        auto const tets =
            subdivide_mesh_for_common_case_3(tetrahedron, {v3, v4, v1, v2}, {ce1}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_3_3)
    {
        auto const& ce1 = cut_edge_keys[e3];
        auto const& cf1 = cut_face_keys[f3];
        auto const& cf2 = cut_face_keys[f4];
        auto const tets =
            subdivide_mesh_for_common_case_3(tetrahedron, {v1, v4, v2, v3}, {ce1}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_3_4)
    {
        auto const& ce1 = cut_edge_keys[e4];
        auto const& cf1 = cut_face_keys[f1];
        auto const& cf2 = cut_face_keys[f3];
        auto const tets =
            subdivide_mesh_for_common_case_3(tetrahedron, {v1, v2, v3, v4}, {ce1}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_3_5)
    {
        auto const& ce1 = cut_edge_keys[e5];
        auto const& cf1 = cut_face_keys[f2];
        auto const& cf2 = cut_face_keys[f1];
        auto const tets =
            subdivide_mesh_for_common_case_3(tetrahedron, {v2, v3, v1, v4}, {ce1}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_3_6)
    {
        auto const& ce1 = cut_edge_keys[e6];
        auto const& cf1 = cut_face_keys[f3];
        auto const& cf2 = cut_face_keys[f2];
        auto const tets =
            subdivide_mesh_for_common_case_3(tetrahedron, {v3, v1, v2, v4}, {ce1}, {cf1, cf2});
        return tets;
    }

    std::byte constexpr case_4_12{0b00000011};
    std::byte constexpr case_4_13{0b00000101};
    std::byte constexpr case_4_14{0b00001001};
    std::byte constexpr case_4_15{0b00010001};
    std::byte constexpr case_4_23{0b00000110};
    std::byte constexpr case_4_25{0b00010010};
    std::byte constexpr case_4_26{0b00100010};
    std::byte constexpr case_4_34{0b00001100};
    std::byte constexpr case_4_36{0b00100100};
    std::byte constexpr case_4_45{0b00011000};
    std::byte constexpr case_4_46{0b00101000};
    std::byte constexpr case_4_56{0b00110000};

    if (edge_intersection_mask == case_4_12)
    {
        auto const& ce1 = cut_edge_keys[e1];
        auto const& ce2 = cut_edge_keys[e2];
        auto const& cf1 = cut_face_keys[f2];
        auto const& cf2 = cut_face_keys[f1];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v1, v3, v4, v2}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_13)
    {
        auto const& ce1 = cut_edge_keys[e3];
        auto const& ce2 = cut_edge_keys[e1];
        auto const& cf1 = cut_face_keys[f1];
        auto const& cf2 = cut_face_keys[f3];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v3, v2, v4, v1}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_14)
    {
        auto const& ce1 = cut_edge_keys[e1];
        auto const& ce2 = cut_edge_keys[e4];
        auto const& cf1 = cut_face_keys[f3];
        auto const& cf2 = cut_face_keys[f4];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v2, v4, v3, v1}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_15)
    {
        auto const& ce1 = cut_edge_keys[e5];
        auto const& ce2 = cut_edge_keys[e1];
        auto const& cf1 = cut_face_keys[f4];
        auto const& cf2 = cut_face_keys[f2];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v4, v1, v3, v2}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_23)
    {
        auto const& ce1 = cut_edge_keys[e2];
        auto const& ce2 = cut_edge_keys[e3];
        auto const& cf1 = cut_face_keys[f3];
        auto const& cf2 = cut_face_keys[f2];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v2, v1, v4, v3}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_25)
    {
        auto const& ce1 = cut_edge_keys[e2];
        auto const& ce2 = cut_edge_keys[e5];
        auto const& cf1 = cut_face_keys[f1];
        auto const& cf2 = cut_face_keys[f4];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v3, v4, v1, v2}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_26)
    {
        auto const& ce1 = cut_edge_keys[e6];
        auto const& ce2 = cut_edge_keys[e2];
        auto const& cf1 = cut_face_keys[f4];
        auto const& cf2 = cut_face_keys[f3];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v4, v2, v1, v3}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_34)
    {
        auto const& ce1 = cut_edge_keys[e4];
        auto const& ce2 = cut_edge_keys[e3];
        auto const& cf1 = cut_face_keys[f4];
        auto const& cf2 = cut_face_keys[f1];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v4, v3, v2, v1}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_36)
    {
        auto const& ce1 = cut_edge_keys[e3];
        auto const& ce2 = cut_edge_keys[e6];
        auto const& cf1 = cut_face_keys[f2];
        auto const& cf2 = cut_face_keys[f4];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v1, v4, v2, v3}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_45)
    {
        auto const& ce1 = cut_edge_keys[e4];
        auto const& ce2 = cut_edge_keys[e5];
        auto const& cf1 = cut_face_keys[f2];
        auto const& cf2 = cut_face_keys[f3];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v1, v2, v3, v4}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_46)
    {
        auto const& ce1 = cut_edge_keys[e6];
        auto const& ce2 = cut_edge_keys[e4];
        auto const& cf1 = cut_face_keys[f1];
        auto const& cf2 = cut_face_keys[f2];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v3, v1, v2, v4}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_56)
    {
        auto const& ce1 = cut_edge_keys[e5];
        auto const& ce2 = cut_edge_keys[e6];
        auto const& cf1 = cut_face_keys[f3];
        auto const& cf2 = cut_face_keys[f1];
        auto const tets =
            subdivide_mesh_for_common_case_4(tetrahedron, {v2, v3, v1, v4}, {ce1, ce2}, {cf1, cf2});
        return tets;
    }

    std::byte constexpr case_5_124{0b00001011};
    std::byte constexpr case_5_126{0b00100011};
    std::byte constexpr case_5_135{0b00010101};
    std::byte constexpr case_5_136{0b00100101};
    std::byte constexpr case_5_146{0b00101001};
    std::byte constexpr case_5_156{0b00110001};
    std::byte constexpr case_5_234{0b00001110};
    std::byte constexpr case_5_235{0b00010110};
    std::byte constexpr case_5_245{0b00011010};
    std::byte constexpr case_5_246{0b00101010};
    std::byte constexpr case_5_345{0b00011100};
    std::byte constexpr case_5_356{0b00110100};

    bool constexpr symmetry{true};
    if (edge_intersection_mask == case_5_124)
    {
        auto const& ce1 = cut_edge_keys[e4];
        auto const& ce2 = cut_edge_keys[e1];
        auto const& ce3 = cut_edge_keys[e2];
        auto const& cf1 = cut_face_keys[f2];
        auto const& cf2 = cut_face_keys[f3];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v4, v2, v3, v1},
            {ce1, ce2, ce3},
            {cf1, cf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_126)
    {
        auto const& ce1 = cut_edge_keys[e1];
        auto const& ce2 = cut_edge_keys[e2];
        auto const& ce3 = cut_edge_keys[e6];
        auto const& cf1 = cut_face_keys[f3];
        auto const& cf2 = cut_face_keys[f1];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v1, v3, v4, v2},
            {ce1, ce2, ce3},
            {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_5_135)
    {
        auto const& ce1 = cut_edge_keys[e5];
        auto const& ce2 = cut_edge_keys[e1];
        auto const& ce3 = cut_edge_keys[e3];
        auto const& cf1 = cut_face_keys[f3];
        auto const& cf2 = cut_face_keys[f2];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v4, v1, v3, v2},
            {ce1, ce2, ce3},
            {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_5_136)
    {
        auto const& ce1 = cut_edge_keys[e6];
        auto const& ce2 = cut_edge_keys[e3];
        auto const& ce3 = cut_edge_keys[e1];
        auto const& cf1 = cut_face_keys[f1];
        auto const& cf2 = cut_face_keys[f2];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v4, v1, v2, v3},
            {ce1, ce2, ce3},
            {cf1, cf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_146)
    {
        auto const& ce1 = cut_edge_keys[e1];
        auto const& ce2 = cut_edge_keys[e4];
        auto const& ce3 = cut_edge_keys[e6];
        auto const& cf1 = cut_face_keys[f2];
        auto const& cf2 = cut_face_keys[f4];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v2, v4, v3, v1},
            {ce1, ce2, ce3},
            {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_5_156)
    {
        auto const& ce1 = cut_edge_keys[e6];
        auto const& ce2 = cut_edge_keys[e5];
        auto const& ce3 = cut_edge_keys[e1];
        auto const& cf1 = cut_face_keys[f4];
        auto const& cf2 = cut_face_keys[f3];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v3, v2, v1, v4},
            {ce1, ce2, ce3},
            {cf1, cf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_234)
    {
        auto const& ce1 = cut_edge_keys[e2];
        auto const& ce2 = cut_edge_keys[e3];
        auto const& ce3 = cut_edge_keys[e4];
        auto const& cf1 = cut_face_keys[f1];
        auto const& cf2 = cut_face_keys[f2];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v2, v1, v4, v3},
            {ce1, ce2, ce3},
            {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_5_235)
    {
        auto const& ce1 = cut_edge_keys[e3];
        auto const& ce2 = cut_edge_keys[e2];
        auto const& ce3 = cut_edge_keys[e5];
        auto const& cf1 = cut_face_keys[f1];
        auto const& cf2 = cut_face_keys[f3];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v1, v2, v4, v3},
            {ce1, ce2, ce3},
            {cf1, cf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_245)
    {
        auto const& ce1 = cut_edge_keys[e4];
        auto const& ce2 = cut_edge_keys[e5];
        auto const& ce3 = cut_edge_keys[e2];
        auto const& cf1 = cut_face_keys[f4];
        auto const& cf2 = cut_face_keys[f3];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v1, v2, v3, v4},
            {ce1, ce2, ce3},
            {cf1, cf2});
        return tets;
    }
    if (edge_intersection_mask == case_5_246)
    {
        auto const& ce1 = cut_edge_keys[e4];
        auto const& ce2 = cut_edge_keys[e6];
        auto const& ce3 = cut_edge_keys[e2];
        auto const& cf1 = cut_face_keys[f4];
        auto const& cf2 = cut_face_keys[f1];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v1, v3, v2, v4},
            {ce1, ce2, ce3},
            {cf1, cf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_345)
    {
        auto const& ce1 = cut_edge_keys[e5];
        auto const& ce2 = cut_edge_keys[e4];
        auto const& ce3 = cut_edge_keys[e3];
        auto const& cf1 = cut_face_keys[f4];
        auto const& cf2 = cut_face_keys[f2];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v2, v1, v3, v4},
            {ce1, ce2, ce3},
            {cf1, cf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_356)
    {
        auto const& ce1 = cut_edge_keys[e5];
        auto const& ce2 = cut_edge_keys[e6];
        auto const& ce3 = cut_edge_keys[e3];
        auto const& cf1 = cut_face_keys[f4];
        auto const& cf2 = cut_face_keys[f1];
        auto const tets = subdivide_mesh_for_common_case_5(
            tetrahedron,
            {v2, v3, v1, v4},
            {ce1, ce2, ce3},
            {cf1, cf2});
        return tets;
    }

    return {};
}

tetrahedral_mesh_cutter_t::subdivided_element_type
tetrahedral_mesh_cutter_t::subdivide_mesh_for_common_case_1(
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<edge_facet_key_type, 3u> const& edge_intersections)
{
    std::uint32_t const _v1 = mesh_.elements().col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = mesh_.elements().col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = mesh_.elements().col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = mesh_.elements().col(tetrahedron)(vertex_ordering[3]);

    auto const ce1_v1 = _v1;
    auto const ce1_v2 = _v4;
    auto const ce2_v1 = _v2;
    auto const ce2_v2 = _v4;
    auto const ce3_v1 = _v3;
    auto const ce3_v2 = _v4;

    auto const& ce1 = cut_edges_[{ce1_v1, ce1_v2}];
    auto const& ce2 = cut_edges_[{ce2_v1, ce2_v2}];
    auto const& ce3 = cut_edges_[{ce3_v1, ce3_v2}];

    std::uint32_t const v1 = _v1;
    std::uint32_t const v2 = _v2;
    std::uint32_t const v3 = _v3;

    bool const should_not_reverse_ce1 = (ce1.v1 == ce1_v1) && (ce1.v2 == ce1_v2);
    bool const should_not_reverse_ce2 = (ce2.v1 == ce2_v1) && (ce2.v2 == ce2_v2);
    bool const should_not_reverse_ce3 = (ce3.v1 == ce3_v1) && (ce3.v2 == ce3_v2);

    std::uint32_t const v5 = should_not_reverse_ce1 ? ce1.vi : ce1.vi + 1u;
    std::uint32_t const v6 = should_not_reverse_ce2 ? ce2.vi : ce2.vi + 1u;
    std::uint32_t const v7 = should_not_reverse_ce3 ? ce3.vi : ce3.vi + 1u;

    std::uint32_t const v4  = _v4;
    std::uint32_t const v5p = should_not_reverse_ce1 ? ce1.vi + 1u : ce1.vi;
    std::uint32_t const v6p = should_not_reverse_ce2 ? ce2.vi + 1u : ce2.vi;
    std::uint32_t const v7p = should_not_reverse_ce3 ? ce3.vi + 1u : ce3.vi;

    std::size_t constexpr num_new_tetrahedra = 4u;
    tetrahedra_type T(4u, num_new_tetrahedra);

    T.col(0u) = tetrahedron_type{v1, v2, v3, v7};
    T.col(1u) = tetrahedron_type{v1, v5, v2, v7};
    T.col(2u) = tetrahedron_type{v6, v5, v7, v2};
    T.col(3u) = tetrahedron_type{v5p, v6p, v7p, v4};

    return T;
}

tetrahedral_mesh_cutter_t::subdivided_element_type
tetrahedral_mesh_cutter_t::subdivide_mesh_for_common_case_2(
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<edge_facet_key_type, 4u> const& edge_intersections)
{
    std::uint32_t const _v1 = mesh_.elements().col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = mesh_.elements().col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = mesh_.elements().col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = mesh_.elements().col(tetrahedron)(vertex_ordering[3]);

    auto const ce1_v1 = _v1;
    auto const ce1_v2 = _v4;
    auto const ce2_v1 = _v3;
    auto const ce2_v2 = _v1;
    auto const ce3_v1 = _v2;
    auto const ce3_v2 = _v3;
    auto const ce4_v1 = _v2;
    auto const ce4_v2 = _v4;

    auto const& ce1 = cut_edges_[{ce1_v1, ce1_v2}];
    auto const& ce2 = cut_edges_[{ce2_v1, ce2_v2}];
    auto const& ce3 = cut_edges_[{ce3_v1, ce3_v2}];
    auto const& ce4 = cut_edges_[{ce4_v1, ce4_v2}];

    bool const should_not_reverse_ce1 = (ce1.v1 == ce1_v1) && (ce1.v2 == ce1_v2);
    bool const should_not_reverse_ce2 = (ce2.v1 == ce2_v2) && (ce2.v2 == ce2_v1);
    bool const should_not_reverse_ce3 = (ce3.v1 == ce3_v1) && (ce3.v2 == ce3_v2);
    bool const should_not_reverse_ce4 = (ce4.v1 == ce4_v1) && (ce4.v2 == ce4_v2);

    std::uint32_t const v3 = _v3;
    std::uint32_t const v4 = _v4;
    std::uint32_t const v5 = should_not_reverse_ce1 ? ce1.vi + 1u : ce1.vi;
    std::uint32_t const v6 = should_not_reverse_ce2 ? ce2.vi + 1u : ce2.vi;
    std::uint32_t const v7 = should_not_reverse_ce3 ? ce3.vi + 1u : ce3.vi;
    std::uint32_t const v8 = should_not_reverse_ce4 ? ce4.vi + 1u : ce4.vi;

    std::uint32_t const v1  = _v1;
    std::uint32_t const v2  = _v2;
    std::uint32_t const v5p = should_not_reverse_ce1 ? ce1.vi : ce1.vi + 1u;
    std::uint32_t const v6p = should_not_reverse_ce2 ? ce2.vi : ce2.vi + 1u;
    std::uint32_t const v7p = should_not_reverse_ce3 ? ce3.vi : ce3.vi + 1u;
    std::uint32_t const v8p = should_not_reverse_ce4 ? ce4.vi : ce4.vi + 1u;

    std::uint32_t constexpr num_new_tetrahedra = 6u;
    tetrahedra_type T(4u, num_new_tetrahedra);

    T.col(0u) = tetrahedron_type{v5, v8, v6, v4};
    T.col(1u) = tetrahedron_type{v6, v8, v7, v4};
    T.col(2u) = tetrahedron_type{v6, v7, v3, v4};

    T.col(3u) = tetrahedron_type{v5p, v6p, v8p, v1};
    T.col(4u) = tetrahedron_type{v1, v2, v6p, v8p};
    T.col(5u) = tetrahedron_type{v2, v7p, v6p, v8p};

    return T;
}

tetrahedral_mesh_cutter_t::subdivided_element_type
tetrahedral_mesh_cutter_t::subdivide_mesh_for_common_case_3(
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<edge_facet_key_type, 1u> const& edge_intersections,
    std::array<triangle_facet_key_type, 2u> const& face_intersections)
{
    std::uint32_t const _v1 = mesh_.elements().col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = mesh_.elements().col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = mesh_.elements().col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = mesh_.elements().col(tetrahedron)(vertex_ordering[3]);

    auto const ce1_v1 = _v1;
    auto const ce1_v2 = _v4;

    auto const& ce1 = cut_edges_[{ce1_v1, ce1_v2}];
    auto const& cf1 = cut_faces_[face_intersections[0]];
    auto const& cf2 = cut_faces_[face_intersections[1]];

    bool const should_not_reverse_ce1 = (ce1.v1 == ce1_v1) && (ce1.v2 == ce1_v2);

    std::uint32_t const v1  = _v1;
    std::uint32_t const v2  = _v2;
    std::uint32_t const v3  = _v3;
    std::uint32_t const v4  = _v4;
    std::uint32_t const v5  = should_not_reverse_ce1 ? ce1.vi : ce1.vi + 1u;
    std::uint32_t const v6  = cf1.vi;
    std::uint32_t const v7  = cf2.vi;
    std::uint32_t const v5p = should_not_reverse_ce1 ? ce1.vi + 1u : ce1.vi;

    std::uint32_t constexpr num_new_tetrahedra = 6u;
    tetrahedra_type T(4u, num_new_tetrahedra);

    T.col(0u) = tetrahedron_type{v5p, v6, v7, v4};
    T.col(1u) = tetrahedron_type{v5, v7, v6, v1};
    T.col(2u) = tetrahedron_type{v6, v3, v7, v4};
    T.col(3u) = tetrahedron_type{v2, v3, v6, v4};
    T.col(4u) = tetrahedron_type{v1, v7, v6, v3};
    T.col(5u) = tetrahedron_type{v1, v2, v3, v6};

    return T;
}

tetrahedral_mesh_cutter_t::subdivided_element_type
tetrahedral_mesh_cutter_t::subdivide_mesh_for_common_case_4(
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<edge_facet_key_type, 2u> const& edge_intersections,
    std::array<triangle_facet_key_type, 2u> const& face_intersections)
{
    std::uint32_t const _v1 = mesh_.elements().col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = mesh_.elements().col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = mesh_.elements().col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = mesh_.elements().col(tetrahedron)(vertex_ordering[3]);

    auto const ce1_v1 = _v1;
    auto const ce1_v2 = _v4;
    auto const ce2_v1 = _v2;
    auto const ce2_v2 = _v4;

    auto const& ce1 = cut_edges_[{ce1_v1, ce1_v2}];
    auto const& ce2 = cut_edges_[{ce2_v1, ce2_v2}];
    auto const& cf1 = cut_faces_[face_intersections[0]];
    auto const& cf2 = cut_faces_[face_intersections[1]];

    bool const should_not_reverse_ce1 = (ce1.v1 == ce1_v1) && (ce1.v2 == ce1_v2);
    bool const should_not_reverse_ce2 = (ce2.v1 == ce2_v1) && (ce2.v2 == ce2_v2);

    std::uint32_t const v1  = _v1;
    std::uint32_t const v2  = _v2;
    std::uint32_t const v3  = _v3;
    std::uint32_t const v4  = _v4;
    std::uint32_t const v5  = should_not_reverse_ce1 ? ce1.vi : ce1.vi + 1u;
    std::uint32_t const v6  = should_not_reverse_ce2 ? ce2.vi : ce2.vi + 1u;
    std::uint32_t const v7  = cf1.vi;
    std::uint32_t const v8  = cf2.vi;
    std::uint32_t const v5p = should_not_reverse_ce1 ? ce1.vi + 1u : ce1.vi;
    std::uint32_t const v6p = should_not_reverse_ce2 ? ce2.vi + 1u : ce2.vi;

    std::uint32_t constexpr num_new_tetrahedra = 8u;
    tetrahedra_type T(4u, num_new_tetrahedra);

    T.col(0u) = tetrahedron_type{v5p, v6p, v7, v4};
    T.col(1u) = tetrahedron_type{v5p, v7, v8, v4};
    T.col(2u) = tetrahedron_type{v4, v7, v8, v3};
    T.col(3u) = tetrahedron_type{v7, v6, v5, v2};
    T.col(4u) = tetrahedron_type{v8, v7, v5, v2};
    T.col(5u) = tetrahedron_type{v8, v7, v2, v3};
    T.col(6u) = tetrahedron_type{v5, v8, v2, v1};
    T.col(7u) = tetrahedron_type{v1, v2, v3, v8};

    return T;
}

tetrahedral_mesh_cutter_t::subdivided_element_type
tetrahedral_mesh_cutter_t::subdivide_mesh_for_common_case_5(
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<edge_facet_key_type, 3u> const& edge_intersections,
    std::array<triangle_facet_key_type, 2u> const& face_intersections,
    bool symmetry)
{
    std::uint32_t const _v1 = mesh_.elements().col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = mesh_.elements().col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = mesh_.elements().col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = mesh_.elements().col(tetrahedron)(vertex_ordering[3]);

    auto const ce1_v1 = _v1;
    auto const ce1_v2 = _v4;
    auto const ce2_v1 = _v2;
    auto const ce2_v2 = _v4;
    auto const ce3_v1 = _v2;
    auto const ce3_v2 = _v3;

    auto const& ce1 = cut_edges_[{ce1_v1, ce1_v2}];
    auto const& ce2 = cut_edges_[{ce2_v1, ce2_v2}];
    auto const& ce3 = cut_edges_[{ce3_v1, ce3_v2}];
    auto const& cf1 = cut_faces_[face_intersections[0]];
    auto const& cf2 = cut_faces_[face_intersections[1]];

    bool const should_not_reverse_ce1 = (ce1.v1 == ce1_v1) && (ce1.v2 == ce1_v2);
    bool const should_not_reverse_ce3 = (ce3.v1 == ce3_v1) && (ce3.v2 == ce3_v2);
    bool const should_not_reverse_ce2 = (ce2.v1 == ce2_v1) && (ce2.v2 == ce2_v2);

    std::uint32_t const v1  = _v1;
    std::uint32_t const v2  = _v2;
    std::uint32_t const v3  = _v3;
    std::uint32_t const v4  = _v4;
    std::uint32_t const v5  = should_not_reverse_ce1 ? ce1.vi : ce1.vi + 1u;
    std::uint32_t const v6  = should_not_reverse_ce2 ? ce2.vi : ce2.vi + 1u;
    std::uint32_t const v7  = should_not_reverse_ce3 ? ce3.vi : ce3.vi + 1u;
    std::uint32_t const v8  = cf1.vi;
    std::uint32_t const v9  = cf2.vi;
    std::uint32_t const v5p = should_not_reverse_ce1 ? ce1.vi + 1u : ce1.vi;
    std::uint32_t const v6p = should_not_reverse_ce2 ? ce2.vi + 1u : ce2.vi;
    std::uint32_t const v7p = should_not_reverse_ce3 ? ce3.vi + 1u : ce3.vi;

    std::uint32_t constexpr num_new_tetrahedra = 9u;
    tetrahedra_type T(4u, num_new_tetrahedra);

    if (symmetry)
    {
        T.col(0u) = tetrahedron_type{v9, v5, v6, v1};
        T.col(1u) = tetrahedron_type{v9, v6, v8, v1};
        T.col(2u) = tetrahedron_type{v1, v3, v8, v9};
        T.col(3u) = tetrahedron_type{v8, v6, v7, v1};
        T.col(4u) = tetrahedron_type{v2, v7, v6, v1};
        T.col(5u) = tetrahedron_type{v5p, v9, v6p, v4};
        T.col(6u) = tetrahedron_type{v9, v3, v6p, v4};
        T.col(7u) = tetrahedron_type{v9, v3, v8, v6p};
        T.col(8u) = tetrahedron_type{v8, v3, v7p, v6p};
    }
    else
    {
        T.col(0u) = tetrahedron_type{v9, v6, v5, v1};
        T.col(1u) = tetrahedron_type{v9, v8, v6, v1};
        T.col(2u) = tetrahedron_type{v1, v8, v3, v9};
        T.col(3u) = tetrahedron_type{v8, v7, v6, v1};
        T.col(4u) = tetrahedron_type{v2, v6, v7, v1};
        T.col(5u) = tetrahedron_type{v5p, v6p, v9, v4};
        T.col(6u) = tetrahedron_type{v9, v6p, v3, v4};
        T.col(7u) = tetrahedron_type{v9, v8, v3, v6p};
        T.col(8u) = tetrahedron_type{v8, v7p, v3, v6p};
    }

    return T;
}

} // namespace cutting
} // namespace physics
} // namespace sbs