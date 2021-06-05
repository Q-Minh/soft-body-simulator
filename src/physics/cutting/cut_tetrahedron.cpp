#include "physics/cutting/cut_tetrahedron.h"

#include "common/primitive.h"

#include <Eigen/LU>
#include <iostream>

namespace sbs {
namespace physics {
namespace cutting {

std::optional<tetrahedron_mesh_cutter_t::subdivided_element_type> cut_tetrahedron(
    common::shared_vertex_mesh_t const& mesh,
    std::uint32_t tetrahedron,
    std::byte const edge_intersection_mask,
    std::array<Eigen::Vector3d, 6u> const& edge_intersections,
    std::array<Eigen::Vector3d, 4u> const& face_intersections)
{
    tetrahedron_mesh_cutter_t cutter{};
    auto const tets = cutter.subdivide_mesh(
        edge_intersection_mask,
        mesh,
        tetrahedron,
        edge_intersections,
        face_intersections);

    return tets;
}

std::optional<tetrahedron_mesh_cutter_t::subdivided_element_type> cut_tetrahedron(
    common::shared_vertex_mesh_t const& mesh,
    std::uint32_t tetrahedron,
    common::shared_vertex_surface_mesh_t const& cutting_surface)
{
    auto const v1 = mesh.elements()(0u, tetrahedron);
    auto const v2 = mesh.elements()(1u, tetrahedron);
    auto const v3 = mesh.elements()(2u, tetrahedron);
    auto const v4 = mesh.elements()(3u, tetrahedron);

    auto const& p1 = mesh.positions().col(v1);
    auto const& p2 = mesh.positions().col(v2);
    auto const& p3 = mesh.positions().col(v3);
    auto const& p4 = mesh.positions().col(v4);

    collision::line_segment_t const e1{p1, p2};
    collision::line_segment_t const e2{p2, p3};
    collision::line_segment_t const e3{p3, p1};
    collision::line_segment_t const e4{p1, p4};
    collision::line_segment_t const e5{p2, p4};
    collision::line_segment_t const e6{p3, p4};

    std::byte edge_intersection_mask{0b00000000};

    std::array<Eigen::Vector3d, 6u> edge_intersections{};
    std::array<Eigen::Vector3d, 4u> face_intersections{};

    /**
     * Find tetrahedron's cut edges
     */
    std::size_t const num_cutting_triangles =
        static_cast<std::size_t>(cutting_surface.triangles().cols());

    for (std::size_t f = 0u; f < num_cutting_triangles; ++f)
    {
        auto const v1 = cutting_surface.triangles()(0u, f);
        auto const v2 = cutting_surface.triangles()(1u, f);
        auto const v3 = cutting_surface.triangles()(2u, f);

        collision::triangle_t const cutting_triangle{
            cutting_surface.vertices().col(v1),
            cutting_surface.vertices().col(v2),
            cutting_surface.vertices().col(v3)};

        auto const e1_intersection = collision::intersect_twoway(e1, cutting_triangle);
        auto const e2_intersection = collision::intersect_twoway(e2, cutting_triangle);
        auto const e3_intersection = collision::intersect_twoway(e3, cutting_triangle);
        auto const e4_intersection = collision::intersect_twoway(e4, cutting_triangle);
        auto const e5_intersection = collision::intersect_twoway(e5, cutting_triangle);
        auto const e6_intersection = collision::intersect_twoway(e6, cutting_triangle);

        if (e1_intersection.has_value())
        {
            edge_intersections[0] = e1_intersection.value();
            edge_intersection_mask |= std::byte{0b00000001};
        }
        if (e2_intersection.has_value())
        {
            edge_intersections[1] = e2_intersection.value();
            edge_intersection_mask |= std::byte{0b00000010};
        }
        if (e3_intersection.has_value())
        {
            edge_intersections[2] = e3_intersection.value();
            edge_intersection_mask |= std::byte{0b00000100};
        }
        if (e4_intersection.has_value())
        {
            edge_intersections[3] = e4_intersection.value();
            edge_intersection_mask |= std::byte{0b00001000};
        }
        if (e5_intersection.has_value())
        {
            edge_intersections[4] = e5_intersection.value();
            edge_intersection_mask |= std::byte{0b00010000};
        }
        if (e6_intersection.has_value())
        {
            edge_intersections[5] = e6_intersection.value();
            edge_intersection_mask |= std::byte{0b00100000};
        }
    }

    collision::triangle_t const f1{p1, p2, p4};
    collision::triangle_t const f2{p2, p3, p4};
    collision::triangle_t const f3{p3, p1, p4};
    collision::triangle_t const f4{p1, p3, p2};

    auto const boundary_edges = cutting_surface.boundary_edges();
    for (auto const& [v1, v2] : boundary_edges)
    {
        collision::line_segment_t const boundary_edge{
            cutting_surface.vertices().col(v1),
            cutting_surface.vertices().col(v2)};

        auto const f1_intersection = collision::intersect_twoway(boundary_edge, f1);
        auto const f2_intersection = collision::intersect_twoway(boundary_edge, f2);
        auto const f3_intersection = collision::intersect_twoway(boundary_edge, f3);
        auto const f4_intersection = collision::intersect_twoway(boundary_edge, f4);

        if (f1_intersection.has_value())
            face_intersections[0] = f1_intersection.value();
        if (f2_intersection.has_value())
            face_intersections[1] = f2_intersection.value();
        if (f3_intersection.has_value())
            face_intersections[2] = f3_intersection.value();
        if (f4_intersection.has_value())
            face_intersections[3] = f4_intersection.value();
    }

    return cut_tetrahedron(
        mesh,
        tetrahedron,
        edge_intersection_mask,
        edge_intersections,
        face_intersections);
}

bool cut_tetrahedral_mesh(
    common::shared_vertex_mesh_t& mesh,
    common::shared_vertex_surface_mesh_t const& cutting_surface)
{
    bool has_mesh_been_cut{false};

    std::uint32_t const previous_number_of_elements =
        static_cast<std::uint32_t>(mesh.elements().cols());
    std::uint32_t const previous_number_of_vertices =
        static_cast<std::uint32_t>(mesh.positions().cols());

    std::uint32_t current_number_of_vertices = previous_number_of_vertices;
    std::uint32_t current_number_of_elements = previous_number_of_elements;

    /**
     * Brute-force search for intersections between the cutting surface and the mesh.
     * Pruning the search space using efficient spatial acceleration structures
     * could be beneficial for performance.
     */
    std::list<std::pair<
        std::uint32_t /* index of tetrahedron that was cut */,
        tetrahedron_mesh_cutter_t::
            subdivided_element_type /* the subdivided tets originating from the cut tetrahedron*/>>
        new_tetrahedra{};

    for (std::uint32_t e = 0u; e < previous_number_of_elements; ++e)
    {
        auto subdivided_tets = cut_tetrahedron(mesh, e, cutting_surface);

        if (!subdivided_tets.has_value())
            continue;

        has_mesh_been_cut = true;

        auto& [P, T, M, V, F] = subdivided_tets.value();

        std::cout << "positions:\n" << P << "\n";
        std::cout << "masses:\n" << M << "\n";
        std::cout << "velocities:\n" << V << "\n";
        std::cout << "forces:\n" << F << "\n";

        std::uint32_t const num_subdivided_tets = static_cast<std::uint32_t>(T.cols());
        /**
         * Check for any vertex index that is >= previous number of elements.
         * Those vertices are newly created vertices. Since we only affect the
         * mesh after having found every cut tetrahedron, every new vertex index
         * must be adjusted to account for all other new vertex indices. For example,
         * for a mesh of size 5 initially, if we cut one of its tets and generate 3 new
         * vertices, then the new vertices will have indices 5,6,7. But if in the
         * same cut, another tet was cut and generated 3 new vertices, we can't
         * also give those the indices 5,6,7. We will have to start from 8 to account
         * for the previously cut tet. So the 6 new vertices will have indices 5,6,7,8,9,10.
         */
        for (std::uint32_t t = 0u; t < num_subdivided_tets; ++t)
        {
            for (std::uint32_t i = 0u; i < 4u; ++i)
            {
                auto const vi = T(i, t);

                if (vi < previous_number_of_vertices)
                    continue;

                auto const offset = vi - previous_number_of_vertices;
                auto const vip    = current_number_of_vertices + offset;
                T(i, t)           = vip;
            }
        }

        std::uint32_t const num_new_vertices = static_cast<std::uint32_t>(P.cols());
        current_number_of_vertices += num_new_vertices;

        /**
         * Here, we subtract the number of subdivided tets by 1, because the initial tetrahedron
         * from which these subdivided tets come from will be removed from the mesh. Thus,
         * we add all subdivided tets to the mesh, but we remove the initial tetrahedron from the
         * mesh.
         */
        current_number_of_elements += num_subdivided_tets - 1u;

        new_tetrahedra.push_back(std::make_pair(e, subdivided_tets.value()));
    }

    mesh.elements().conservativeResize(4u, current_number_of_elements);
    mesh.positions().conservativeResize(3u, current_number_of_vertices);
    mesh.masses().conservativeResize(current_number_of_vertices);
    mesh.velocities().conservativeResize(3u, current_number_of_vertices);
    mesh.forces().conservativeResize(3u, current_number_of_vertices);

    int e_offset = static_cast<int>(previous_number_of_elements);
    int v_offset = static_cast<int>(previous_number_of_vertices);

    for (auto const& [cut_tetrahedron, subdivided_tets] : new_tetrahedra)
    {
        auto const& [P, T, M, V, F] = subdivided_tets;

        std::cout << "positions:\n" << P << "\n";
        std::cout << "masses:\n" << M << "\n";
        std::cout << "velocities:\n" << V << "\n";
        std::cout << "forces:\n" << F << "\n";

        mesh.positions().block(0, v_offset, P.rows(), P.cols())  = P;
        mesh.masses().block(0, v_offset, M.rows(), M.cols())     = M;
        mesh.velocities().block(0, v_offset, V.rows(), V.cols()) = V;
        mesh.forces().block(0, v_offset, F.rows(), F.cols())     = F;

        v_offset += P.cols();

        int const num_tetrahedra_to_append = T.cols() - 1;

        /**
         * We replace the cut tetrahedron by the first tetrahedron of the
         * subdivided tet.
         */
        mesh.elements().col(cut_tetrahedron) = T.col(0);
        /**
         * We then append the second to last subdivided tets the
         * our mesh elements.
         */
        mesh.elements().block(0, e_offset, T.rows(), num_tetrahedra_to_append) =
            T.block(0, 1, T.rows(), num_tetrahedra_to_append);

        e_offset += num_tetrahedra_to_append;
    }

    return has_mesh_been_cut;
}

std::optional<tetrahedron_mesh_cutter_t::subdivided_element_type>
tetrahedron_mesh_cutter_t::subdivide_mesh(
    std::byte const& edge_intersection_mask,
    common::shared_vertex_mesh_t const& mesh,
    std::uint32_t tetrahedron,
    std::array<Eigen::Vector3d, 6u> const& edge_intersection_points,
    std::array<Eigen::Vector3d, 4u> const& face_intersection_points)
{
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
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pe2 = edge_intersection_points[e2];
        auto const& pe3 = edge_intersection_points[e5];
        auto const tets =
            subdivide_mesh_for_common_case_1(mesh, tetrahedron, {v1, v3, v4, v2}, {pe1, pe2, pe3});
        return tets;
    }
    if (edge_intersection_mask == case_1_134)
    {
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pe2 = edge_intersection_points[e4];
        auto const& pe3 = edge_intersection_points[e3];
        auto const tets =
            subdivide_mesh_for_common_case_1(mesh, tetrahedron, {v2, v4, v3, v1}, {pe1, pe2, pe3});
        return tets;
    }
    if (edge_intersection_mask == case_1_236)
    {
        auto const& pe1 = edge_intersection_points[e2];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pe3 = edge_intersection_points[e6];
        auto const tets =
            subdivide_mesh_for_common_case_1(mesh, tetrahedron, {v2, v1, v4, v3}, {pe1, pe2, pe3});
        return tets;
    }
    if (edge_intersection_mask == case_1_456)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e5];
        auto const& pe3 = edge_intersection_points[e6];
        auto const tets =
            subdivide_mesh_for_common_case_1(mesh, tetrahedron, {v1, v2, v3, v4}, {pe1, pe2, pe3});
        return tets;
    }

    std::byte constexpr case_2_1246{0b00101011};
    std::byte constexpr case_2_1356{0b00110101};
    std::byte constexpr case_2_2345{0b00011110};

    if (edge_intersection_mask == case_2_1246)
    {
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pe2 = edge_intersection_points[e2];
        auto const& pe3 = edge_intersection_points[e6];
        auto const& pe4 = edge_intersection_points[e4];
        auto const tets = subdivide_mesh_for_common_case_2(
            mesh,
            tetrahedron,
            {v2, v4, v3, v1},
            {pe1, pe2, pe3, pe4});
        return tets;
    }
    if (edge_intersection_mask == case_2_1356)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e1];
        auto const& pe3 = edge_intersection_points[e3];
        auto const& pe4 = edge_intersection_points[e6];
        auto const tets = subdivide_mesh_for_common_case_2(
            mesh,
            tetrahedron,
            {v2, v3, v1, v4},
            {pe1, pe2, pe3, pe4});
        return tets;
    }
    if (edge_intersection_mask == case_2_2345)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pe3 = edge_intersection_points[e2];
        auto const& pe4 = edge_intersection_points[e5];
        auto const tets = subdivide_mesh_for_common_case_2(
            mesh,
            tetrahedron,
            {v1, v2, v3, v4},
            {pe1, pe2, pe3, pe4});
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
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f1];
        auto const tets = subdivide_mesh_for_common_case_3(
            mesh,
            tetrahedron,
            {v1, v3, v4, v2},
            {pe1},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_3_2)
    {
        auto const& pe1 = edge_intersection_points[e2];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f4];
        auto const tets = subdivide_mesh_for_common_case_3(
            mesh,
            tetrahedron,
            {v3, v4, v1, v2},
            {pe1},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_3_3)
    {
        auto const& pe1 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f4];
        auto const tets = subdivide_mesh_for_common_case_3(
            mesh,
            tetrahedron,
            {v1, v4, v2, v3},
            {pe1},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_3_4)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f3];
        auto const tets = subdivide_mesh_for_common_case_3(
            mesh,
            tetrahedron,
            {v1, v2, v3, v4},
            {pe1},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_3_5)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f1];
        auto const tets = subdivide_mesh_for_common_case_3(
            mesh,
            tetrahedron,
            {v2, v3, v1, v4},
            {pe1},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_3_6)
    {
        auto const& pe1 = edge_intersection_points[e6];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f2];
        auto const tets = subdivide_mesh_for_common_case_3(
            mesh,
            tetrahedron,
            {v3, v1, v2, v4},
            {pe1},
            {pf1, pf2});
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
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pe2 = edge_intersection_points[e2];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f1];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v1, v3, v4, v2},
            {pe1, pe2},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_13)
    {
        auto const& pe1 = edge_intersection_points[e3];
        auto const& pe2 = edge_intersection_points[e1];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f3];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v3, v2, v4, v1},
            {pe1, pe2},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_14)
    {
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pe2 = edge_intersection_points[e4];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f4];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v2, v4, v3, v1},
            {pe1, pe2},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_15)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e1];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f2];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v4, v1, v3, v2},
            {pe1, pe2},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_23)
    {
        auto const& pe1 = edge_intersection_points[e2];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f2];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v2, v1, v4, v3},
            {pe1, pe2},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_25)
    {
        auto const& pe1 = edge_intersection_points[e2];
        auto const& pe2 = edge_intersection_points[e5];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f4];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v3, v4, v1, v2},
            {pe1, pe2},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_26)
    {
        auto const& pe1 = edge_intersection_points[e6];
        auto const& pe2 = edge_intersection_points[e2];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f3];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v4, v2, v1, v3},
            {pe1, pe2},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_34)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f1];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v4, v3, v2, v1},
            {pe1, pe2},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_36)
    {
        auto const& pe1 = edge_intersection_points[e3];
        auto const& pe2 = edge_intersection_points[e6];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f4];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v1, v4, v2, v3},
            {pe1, pe2},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_45)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e5];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f3];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v1, v2, v3, v4},
            {pe1, pe2},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_46)
    {
        auto const& pe1 = edge_intersection_points[e6];
        auto const& pe2 = edge_intersection_points[e4];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f2];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v3, v1, v2, v4},
            {pe1, pe2},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_4_56)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e6];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f1];
        auto const tets = subdivide_mesh_for_common_case_4(
            mesh,
            tetrahedron,
            {v2, v3, v1, v4},
            {pe1, pe2},
            {pf1, pf2});
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
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e1];
        auto const& pe3 = edge_intersection_points[e2];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f3];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v4, v2, v3, v1},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_126)
    {
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pe2 = edge_intersection_points[e2];
        auto const& pe3 = edge_intersection_points[e6];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f1];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v1, v3, v4, v2},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_5_135)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e1];
        auto const& pe3 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f2];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v4, v1, v3, v2},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_5_136)
    {
        auto const& pe1 = edge_intersection_points[e6];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pe3 = edge_intersection_points[e1];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f2];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v4, v1, v2, v3},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_146)
    {
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pe2 = edge_intersection_points[e4];
        auto const& pe3 = edge_intersection_points[e6];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f4];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v2, v4, v3, v1},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_5_156)
    {
        auto const& pe1 = edge_intersection_points[e6];
        auto const& pe2 = edge_intersection_points[e5];
        auto const& pe3 = edge_intersection_points[e1];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f3];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v3, v2, v1, v4},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_234)
    {
        auto const& pe1 = edge_intersection_points[e2];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pe3 = edge_intersection_points[e4];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f2];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v2, v1, v4, v3},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_5_235)
    {
        auto const& pe1 = edge_intersection_points[e3];
        auto const& pe2 = edge_intersection_points[e2];
        auto const& pe3 = edge_intersection_points[e5];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f3];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v1, v2, v4, v3},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_245)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e5];
        auto const& pe3 = edge_intersection_points[e2];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f3];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v1, v2, v3, v4},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return tets;
    }
    if (edge_intersection_mask == case_5_246)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e6];
        auto const& pe3 = edge_intersection_points[e2];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f1];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v1, v3, v2, v4},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_345)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e4];
        auto const& pe3 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f2];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v2, v1, v3, v4},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return tets;
    }
    if (edge_intersection_mask == case_5_356)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e6];
        auto const& pe3 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f1];
        auto const tets = subdivide_mesh_for_common_case_5(
            mesh,
            tetrahedron,
            {v2, v3, v1, v4},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return tets;
    }

    return {};
}

tetrahedron_mesh_cutter_t::subdivided_element_type
tetrahedron_mesh_cutter_t::subdivide_mesh_for_common_case_1(
    common::shared_vertex_mesh_t const& mesh,
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<Eigen::Vector3d, 3u> const& edge_intersection_points)
{
    std::uint32_t const num_vertices = static_cast<std::uint32_t>(mesh.positions().cols());

    std::uint32_t const _v1 = mesh.elements().col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = mesh.elements().col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = mesh.elements().col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = mesh.elements().col(tetrahedron)(vertex_ordering[3]);

    std::uint32_t const v1 = _v1;
    std::uint32_t const v2 = _v2;
    std::uint32_t const v3 = _v3;

    std::uint32_t constexpr _v5  = 0u;
    std::uint32_t constexpr _v6  = 1u;
    std::uint32_t constexpr _v7  = 2u;
    std::uint32_t constexpr _v5p = 3u;
    std::uint32_t constexpr _v6p = 4u;
    std::uint32_t constexpr _v7p = 5u;

    std::uint32_t const v5 = num_vertices + _v5;
    std::uint32_t const v6 = num_vertices + _v6;
    std::uint32_t const v7 = num_vertices + _v7;

    std::uint32_t const v4  = _v4;
    std::uint32_t const v5p = num_vertices + _v5p;
    std::uint32_t const v6p = num_vertices + _v6p;
    std::uint32_t const v7p = num_vertices + _v7p;

    std::uint32_t constexpr num_new_vertices = 6u;
    positions_type P(3u, num_new_vertices);
    masses_type M(num_new_vertices);
    velocities_type V(3u, num_new_vertices);
    forces_type F(3u, num_new_vertices);

    P.col(_v5) = edge_intersection_points[0];
    P.col(_v6) = edge_intersection_points[1];
    P.col(_v7) = edge_intersection_points[2];

    P.col(_v5p) = P.col(_v5);
    P.col(_v6p) = P.col(_v6);
    P.col(_v7p) = P.col(_v7);

    /**
     * Cache interpolation coefficients
     */
    double const t1 = (edge_intersection_points[0].x() - mesh.positions().col(_v1).x()) /
                      (mesh.positions().col(_v4).x() - mesh.positions().col(_v1).x());
    double const t2 = (edge_intersection_points[1].x() - mesh.positions().col(_v2).x()) /
                      (mesh.positions().col(_v4).x() - mesh.positions().col(_v2).x());
    double const t3 = (edge_intersection_points[2].x() - mesh.positions().col(_v3).x()) /
                      (mesh.positions().col(_v4).x() - mesh.positions().col(_v3).x());

    /**
     * Transfer forces to subdivided mesh elements
     */
    F.col(_v5) = (1 - t1) * mesh.forces().col(_v1) + t1 * mesh.forces().col(_v4);
    F.col(_v6) = (1 - t2) * mesh.forces().col(_v2) + t2 * mesh.forces().col(_v4);
    F.col(_v7) = (1 - t3) * mesh.forces().col(_v3) + t3 * mesh.forces().col(_v4);

    F.col(_v5p) = F.col(_v5);
    F.col(_v6p) = F.col(_v6);
    F.col(_v7p) = F.col(_v7);

    /**
     * Transfer velocities to subdivided mesh elements
     */
    V.col(_v5) = (1 - t1) * mesh.velocities().col(_v1) + t1 * mesh.velocities().col(_v4);
    V.col(_v6) = (1 - t2) * mesh.velocities().col(_v2) + t2 * mesh.velocities().col(_v4);
    V.col(_v7) = (1 - t3) * mesh.velocities().col(_v3) + t3 * mesh.velocities().col(_v4);

    V.col(_v5p) = V.col(_v5);
    V.col(_v6p) = V.col(_v6);
    V.col(_v7p) = V.col(_v7);

    /**
     * Transfer masses to subdivided mesh elements
     */
    M(_v5) = (1 - t1) * mesh.masses()(_v1) + t1 * mesh.masses()(_v4);
    M(_v6) = (1 - t2) * mesh.masses()(_v2) + t2 * mesh.masses()(_v4);
    M(_v7) = (1 - t3) * mesh.masses()(_v3) + t3 * mesh.masses()(_v4);

    M(_v5p) = M(_v5);
    M(_v6p) = M(_v6);
    M(_v7p) = M(_v7);

    std::size_t constexpr num_new_tetrahedra = 4u;
    tetrahedra_type T(4u, num_new_tetrahedra);

    T.col(0u) = tetrahedron_type{v1, v2, v3, v7};
    T.col(1u) = tetrahedron_type{v1, v5, v2, v7};
    T.col(2u) = tetrahedron_type{v6, v5, v7, v2};
    T.col(3u) = tetrahedron_type{v5p, v6p, v7p, v4};

    return std::make_tuple(P, T, M, V, F);
}

tetrahedron_mesh_cutter_t::subdivided_element_type
tetrahedron_mesh_cutter_t::subdivide_mesh_for_common_case_2(
    common::shared_vertex_mesh_t const& mesh,
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<Eigen::Vector3d, 4u> const& edge_intersection_points)
{
    std::uint32_t const num_vertices = static_cast<std::uint32_t>(mesh.positions().cols());

    std::uint32_t const _v1 = mesh.elements().col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = mesh.elements().col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = mesh.elements().col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = mesh.elements().col(tetrahedron)(vertex_ordering[3]);

    std::uint32_t constexpr _v5  = 0u;
    std::uint32_t constexpr _v6  = 1u;
    std::uint32_t constexpr _v7  = 2u;
    std::uint32_t constexpr _v8  = 3u;
    std::uint32_t constexpr _v5p = 4u;
    std::uint32_t constexpr _v6p = 5u;
    std::uint32_t constexpr _v7p = 6u;
    std::uint32_t constexpr _v8p = 7u;

    std::uint32_t const v3 = _v3;
    std::uint32_t const v4 = _v4;
    std::uint32_t const v5 = num_vertices + _v5;
    std::uint32_t const v6 = num_vertices + _v6;
    std::uint32_t const v7 = num_vertices + _v7;
    std::uint32_t const v8 = num_vertices + _v8;

    std::uint32_t const v1  = _v1;
    std::uint32_t const v2  = _v2;
    std::uint32_t const v5p = num_vertices + _v5p;
    std::uint32_t const v6p = num_vertices + _v6p;
    std::uint32_t const v7p = num_vertices + _v7p;
    std::uint32_t const v8p = num_vertices + _v8p;

    std::uint32_t constexpr num_new_vertices = 8u;
    positions_type P(3u, num_new_vertices);
    masses_type M(num_new_vertices);
    velocities_type V(3u, num_new_vertices);
    forces_type F(3u, num_new_vertices);

    P.col(_v5) = edge_intersection_points[0];
    P.col(_v6) = edge_intersection_points[1];
    P.col(_v7) = edge_intersection_points[2];
    P.col(_v8) = edge_intersection_points[3];

    P.col(_v5p) = P.col(_v5);
    P.col(_v6p) = P.col(_v6);
    P.col(_v7p) = P.col(_v7);
    P.col(_v8p) = P.col(_v8);

    /**
     * Cache interpolation coefficients
     */
    double const t1 = (edge_intersection_points[0].x() - mesh.positions().col(_v1).x()) /
                      (mesh.positions().col(_v4).x() - mesh.positions().col(_v1).x());
    double const t2 = (edge_intersection_points[1].x() - mesh.positions().col(_v1).x()) /
                      (mesh.positions().col(_v3).x() - mesh.positions().col(_v1).x());
    double const t3 = (edge_intersection_points[2].x() - mesh.positions().col(_v2).x()) /
                      (mesh.positions().col(_v3).x() - mesh.positions().col(_v2).x());
    double const t4 = (edge_intersection_points[3].x() - mesh.positions().col(_v2).x()) /
                      (mesh.positions().col(_v4).x() - mesh.positions().col(_v2).x());

    /**
     * Transfer forces to subdivided mesh elements
     */
    F.col(_v5) = (1 - t1) * mesh.forces().col(_v1) + t1 * mesh.forces().col(_v4);
    F.col(_v6) = (1 - t2) * mesh.forces().col(_v1) + t2 * mesh.forces().col(_v3);
    F.col(_v7) = (1 - t3) * mesh.forces().col(_v2) + t3 * mesh.forces().col(_v3);
    F.col(_v8) = (1 - t4) * mesh.forces().col(_v2) + t4 * mesh.forces().col(_v4);

    F.col(_v5p) = F.col(_v5);
    F.col(_v6p) = F.col(_v6);
    F.col(_v7p) = F.col(_v7);
    F.col(_v8p) = F.col(_v8);

    /**
     * Transfer velocities to subdivided mesh elements
     */
    V.col(_v5) = (1 - t1) * mesh.velocities().col(_v1) + t1 * mesh.velocities().col(_v4);
    V.col(_v6) = (1 - t2) * mesh.velocities().col(_v1) + t2 * mesh.velocities().col(_v3);
    V.col(_v7) = (1 - t3) * mesh.velocities().col(_v2) + t3 * mesh.velocities().col(_v3);
    V.col(_v8) = (1 - t4) * mesh.velocities().col(_v2) + t4 * mesh.velocities().col(_v4);

    V.col(_v5p) = V.col(_v5);
    V.col(_v6p) = V.col(_v6);
    V.col(_v7p) = V.col(_v7);
    V.col(_v8p) = V.col(_v8);

    /**
     * Transfer masses to subdivided mesh elements
     */
    M(_v5) = (1 - t1) * mesh.masses()(_v1) + t1 * mesh.masses()(_v4);
    M(_v6) = (1 - t2) * mesh.masses()(_v1) + t2 * mesh.masses()(_v3);
    M(_v7) = (1 - t3) * mesh.masses()(_v2) + t3 * mesh.masses()(_v3);
    M(_v8) = (1 - t4) * mesh.masses()(_v2) + t4 * mesh.masses()(_v4);

    M(_v5p) = M(_v5);
    M(_v6p) = M(_v6);
    M(_v7p) = M(_v7);
    M(_v8p) = M(_v8);

    std::uint32_t constexpr num_new_tetrahedra = 6u;
    tetrahedra_type T(4u, num_new_tetrahedra);

    T.col(0u) = tetrahedron_type{v5, v8, v6, v4};
    T.col(1u) = tetrahedron_type{v6, v8, v7, v4};
    T.col(2u) = tetrahedron_type{v6, v7, v3, v4};

    T.col(3u) = tetrahedron_type{v5p, v6p, v8p, v1};
    T.col(4u) = tetrahedron_type{v1, v2, v6p, v8p};
    T.col(5u) = tetrahedron_type{v2, v7p, v6p, v8p};

    return std::make_tuple(P, T, M, V, F);
}

tetrahedron_mesh_cutter_t::subdivided_element_type
tetrahedron_mesh_cutter_t::subdivide_mesh_for_common_case_3(
    common::shared_vertex_mesh_t const& mesh,
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<Eigen::Vector3d, 1u> const& edge_intersection_points,
    std::array<Eigen::Vector3d, 2u> const& face_intersection_points)
{
    std::uint32_t const num_vertices = static_cast<std::uint32_t>(mesh.positions().cols());

    std::uint32_t const _v1 = mesh.elements().col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = mesh.elements().col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = mesh.elements().col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = mesh.elements().col(tetrahedron)(vertex_ordering[3]);

    std::uint32_t constexpr num_new_vertices = 4u;

    std::uint32_t constexpr _v5  = 0u;
    std::uint32_t constexpr _v6  = 1u;
    std::uint32_t constexpr _v7  = 2u;
    std::uint32_t constexpr _v5p = 3u;

    std::uint32_t const v1  = _v1;
    std::uint32_t const v2  = _v2;
    std::uint32_t const v3  = _v3;
    std::uint32_t const v4  = _v4;
    std::uint32_t const v5  = num_vertices + _v5;
    std::uint32_t const v6  = num_vertices + _v6;
    std::uint32_t const v7  = num_vertices + _v7;
    std::uint32_t const v5p = num_vertices + _v5p;

    positions_type P(3u, num_new_vertices);
    masses_type M(num_new_vertices);
    velocities_type V(3u, num_new_vertices);
    forces_type F(3u, num_new_vertices);

    P.col(_v5) = edge_intersection_points[0];
    P.col(_v6) = face_intersection_points[0];
    P.col(_v7) = face_intersection_points[1];

    P.col(_v5p) = P.col(_v5);

    /**
     * Cache interpolation coefficients
     */
    double const t1  = (edge_intersection_points[0].x() - mesh.positions().col(_v1).x()) /
                      (mesh.positions().col(_v4).x() - mesh.positions().col(_v1).x());
    auto const [bu2, bv2, bw2] = common::barycentric_coordinates(
        mesh.positions().col(_v1),
        mesh.positions().col(_v2),
        mesh.positions().col(_v4),
        face_intersection_points[0]);
    auto const [bu3, bv3, bw3] = common::barycentric_coordinates(
        mesh.positions().col(_v1),
        mesh.positions().col(_v4),
        mesh.positions().col(_v3),
        face_intersection_points[1]);

    /**
     * Transfer forces to subdivided mesh elements
     */
    F.col(_v5) = (1 - t1) * mesh.forces().col(_v1) + t1 * mesh.forces().col(_v4);
    F.col(_v6) =
        bu2 * mesh.forces().col(_v1) + bv2 * mesh.forces().col(_v2) + bw2 * mesh.forces().col(_v4);
    F.col(_v7) =
        bu3 * mesh.forces().col(_v1) + bv3 * mesh.forces().col(_v4) + bw3 * mesh.forces().col(_v3);

    F.col(_v5p) = F.col(_v5);

    /**
     * Transfer velocities to subdivided mesh elements
     */
    V.col(_v5) = (1 - t1) * mesh.velocities().col(_v1) + t1 * mesh.velocities().col(_v4);
    V.col(_v6) = bu2 * mesh.velocities().col(_v1) + bv2 * mesh.velocities().col(_v2) +
                 bw2 * mesh.velocities().col(_v4);
    V.col(_v7) = bu3 * mesh.velocities().col(_v1) + bv3 * mesh.velocities().col(_v4) +
                 bw3 * mesh.velocities().col(_v3);

    V.col(_v5p) = V.col(_v5);

    /**
     * Transfer masses to subdivided mesh elements
     */
    M(_v5) = (1 - t1) * mesh.masses()(_v1) + t1 * mesh.masses()(_v4);
    M(_v6) = bu2 * mesh.masses()(_v1) + bv2 * mesh.masses()(_v2) + bw2 * (_v4);
    M(_v7) = bu3 * mesh.masses()(_v1) + bv3 * mesh.masses()(_v4) + bw3 * (_v3);

    M(_v5p) = M(_v5);

    std::uint32_t constexpr num_new_tetrahedra = 6u;
    tetrahedra_type T(4u, num_new_tetrahedra);

    T.col(0u) = tetrahedron_type{v5p, v6, v7, v4};
    T.col(1u) = tetrahedron_type{v5, v7, v6, v1};
    T.col(2u) = tetrahedron_type{v6, v3, v7, v4};
    T.col(3u) = tetrahedron_type{v2, v3, v6, v4};
    T.col(4u) = tetrahedron_type{v1, v7, v6, v3};
    T.col(5u) = tetrahedron_type{v1, v2, v3, v6};

    return std::make_tuple(P, T, M, V, F);
}

tetrahedron_mesh_cutter_t::subdivided_element_type
tetrahedron_mesh_cutter_t::subdivide_mesh_for_common_case_4(
    common::shared_vertex_mesh_t const& mesh,
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<Eigen::Vector3d, 2u> const& edge_intersection_points,
    std::array<Eigen::Vector3d, 2u> const& face_intersection_points)
{
    std::uint32_t const num_vertices = static_cast<std::uint32_t>(mesh.positions().cols());

    std::uint32_t const _v1 = mesh.elements().col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = mesh.elements().col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = mesh.elements().col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = mesh.elements().col(tetrahedron)(vertex_ordering[3]);

    std::uint32_t constexpr num_new_vertices = 6u;

    std::uint32_t constexpr _v5  = 0u;
    std::uint32_t constexpr _v6  = 1u;
    std::uint32_t constexpr _v7  = 2u;
    std::uint32_t constexpr _v8  = 3u;
    std::uint32_t constexpr _v5p = 4u;
    std::uint32_t constexpr _v6p = 5u;

    std::uint32_t const v1  = _v1;
    std::uint32_t const v2  = _v2;
    std::uint32_t const v3  = _v3;
    std::uint32_t const v4  = _v4;
    std::uint32_t const v5  = num_vertices + _v5;
    std::uint32_t const v6  = num_vertices + _v6;
    std::uint32_t const v7  = num_vertices + _v7;
    std::uint32_t const v8  = num_vertices + _v8;
    std::uint32_t const v5p = num_vertices + _v5p;
    std::uint32_t const v6p = num_vertices + _v6p;

    positions_type P(3u, num_new_vertices);
    masses_type M(num_new_vertices);
    velocities_type V(3u, num_new_vertices);
    forces_type F(3u, num_new_vertices);

    P.col(_v5) = edge_intersection_points[0];
    P.col(_v6) = edge_intersection_points[1];
    P.col(_v7) = face_intersection_points[0];
    P.col(_v8) = face_intersection_points[1];

    P.col(_v5p) = P.col(_v5);
    P.col(_v6p) = P.col(_v6);

    /**
     * Cache interpolation coefficients
     */
    double const t1 = (edge_intersection_points[0].x() - mesh.positions().col(_v1).x()) /
                      (mesh.positions().col(_v4).x() - mesh.positions().col(_v1).x());
    double const t2 = (edge_intersection_points[1].x() - mesh.positions().col(_v2).x()) /
                      (mesh.positions().col(_v4).x() - mesh.positions().col(_v2).x());
    auto const [bu3, bv3, bw3] = common::barycentric_coordinates(
        mesh.positions().col(_v2),
        mesh.positions().col(_v3),
        mesh.positions().col(_v4),
        face_intersection_points[0]);
    auto const [bu4, bv4, bw4] = common::barycentric_coordinates(
        mesh.positions().col(_v1),
        mesh.positions().col(_v4),
        mesh.positions().col(_v3),
        face_intersection_points[1]);

    /**
     * Transfer forces to subdivided mesh elements
     */
    F.col(_v5) = (1 - t1) * mesh.forces().col(_v1) + t1 * mesh.forces().col(_v4);
    F.col(_v6) = (1 - t2) * mesh.forces().col(_v2) + t2 * mesh.forces().col(_v4);
    F.col(_v7) =
        bu3 * mesh.forces().col(_v2) + bv3 * mesh.forces().col(_v3) + bw3 * mesh.forces().col(_v4);
    F.col(_v8) =
        bu4 * mesh.forces().col(_v1) + bv4 * mesh.forces().col(_v4) + bw4 * mesh.forces().col(_v3);

    F.col(_v5p) = F.col(_v5);
    F.col(_v6p) = F.col(_v6);

    /**
     * Transfer velocities to subdivided mesh elements
     */
    V.col(_v5) = (1 - t1) * mesh.velocities().col(_v1) + t1 * mesh.velocities().col(_v4);
    V.col(_v6) = (1 - t2) * mesh.velocities().col(_v2) + t2 * mesh.velocities().col(_v4);
    V.col(_v7) = bu3 * mesh.velocities().col(_v2) + bv3 * mesh.velocities().col(_v3) +
                 bw3 * mesh.velocities().col(_v4);
    V.col(_v8) = bu4 * mesh.velocities().col(_v1) + bv4 * mesh.velocities().col(_v4) +
                 bw4 * mesh.velocities().col(_v3);

    V.col(_v5p) = V.col(_v5);
    V.col(_v6p) = V.col(_v6);

    /**
     * Transfer masses to subdivided mesh elements
     */
    M(_v5) = (1 - t1) * mesh.masses()(_v1) + t1 * mesh.masses()(_v4);
    M(_v6) = (1 - t2) * mesh.masses()(_v2) + t2 * mesh.masses()(_v4);
    M(_v7) = bu3 * mesh.masses()(_v2) + bv3 * mesh.masses()(_v3) + bw3 * mesh.masses()(_v4);
    M(_v8) = bu4 * mesh.masses()(_v1) + bv4 * mesh.masses()(_v4) + bw4 * mesh.masses()(_v3);

    M(_v5p) = M(_v5);
    M(_v6p) = M(_v6);

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

    return std::make_tuple(P, T, M, V, F);
}

tetrahedron_mesh_cutter_t::subdivided_element_type
tetrahedron_mesh_cutter_t::subdivide_mesh_for_common_case_5(
    common::shared_vertex_mesh_t const& mesh,
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<Eigen::Vector3d, 3u> const& edge_intersection_points,
    std::array<Eigen::Vector3d, 2u> const& face_intersection_points,
    bool symmetry)
{
    std::uint32_t const num_vertices = static_cast<std::uint32_t>(mesh.positions().cols());

    std::uint32_t const _v1 = mesh.elements().col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = mesh.elements().col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = mesh.elements().col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = mesh.elements().col(tetrahedron)(vertex_ordering[3]);

    std::uint32_t constexpr num_new_vertices = 8u;

    std::uint32_t constexpr _v5  = 0u;
    std::uint32_t constexpr _v6  = 1u;
    std::uint32_t constexpr _v7  = 2u;
    std::uint32_t constexpr _v8  = 3u;
    std::uint32_t constexpr _v9  = 4u;
    std::uint32_t constexpr _v5p = 5u;
    std::uint32_t constexpr _v6p = 6u;
    std::uint32_t constexpr _v7p = 7u;

    std::uint32_t const v1  = _v1;
    std::uint32_t const v2  = _v2;
    std::uint32_t const v3  = _v3;
    std::uint32_t const v4  = _v4;
    std::uint32_t const v5  = num_vertices + _v5;
    std::uint32_t const v6  = num_vertices + _v6;
    std::uint32_t const v7  = num_vertices + _v7;
    std::uint32_t const v8  = num_vertices + _v8;
    std::uint32_t const v9  = num_vertices + _v9;
    std::uint32_t const v5p = num_vertices + _v5p;
    std::uint32_t const v6p = num_vertices + _v6p;
    std::uint32_t const v7p = num_vertices + _v7p;

    positions_type P(3u, num_new_vertices);
    masses_type M(num_new_vertices);
    velocities_type V(3u, num_new_vertices);
    forces_type F(3u, num_new_vertices);

    P.col(_v5) = edge_intersection_points[0];
    P.col(_v6) = edge_intersection_points[1];
    P.col(_v7) = edge_intersection_points[2];
    P.col(_v8) = face_intersection_points[0];
    P.col(_v9) = face_intersection_points[1];

    P.col(_v5p) = P.col(_v5);
    P.col(_v6p) = P.col(_v6);
    P.col(_v7p) = P.col(_v7);

    /**
     * Cache interpolation coefficients
     */
    double const t1 = (edge_intersection_points[0].x() - mesh.positions().col(_v1).x()) /
                      (mesh.positions().col(_v4).x() - mesh.positions().col(_v1).x());
    double const t2 = (edge_intersection_points[1].x() - mesh.positions().col(_v2).x()) /
                      (mesh.positions().col(_v4).x() - mesh.positions().col(_v2).x());
    double const t3 = (edge_intersection_points[2].x() - mesh.positions().col(_v2).x()) /
                      (mesh.positions().col(_v3).x() - mesh.positions().col(_v2).x());
    auto const [bu4, bv4, bw4] = common::barycentric_coordinates(
        mesh.positions().col(_v1),
        mesh.positions().col(_v3),
        mesh.positions().col(_v2),
        face_intersection_points[0]);
    auto const [bu5, bv5, bw5] = common::barycentric_coordinates(
        mesh.positions().col(_v1),
        mesh.positions().col(_v4),
        mesh.positions().col(_v3),
        face_intersection_points[1]);

    /**
     * Transfer forces to subdivided mesh elements
     */
    F.col(_v5) = (1 - t1) * mesh.forces().col(_v1) + t1 * mesh.forces().col(_v4);
    F.col(_v6) = (1 - t2) * mesh.forces().col(_v2) + t2 * mesh.forces().col(_v4);
    F.col(_v7) = (1 - t3) * mesh.forces().col(_v2) + t3 * mesh.forces().col(_v3);
    F.col(_v8) =
        bu4 * mesh.forces().col(_v1) + bv4 * mesh.forces().col(_v3) + bw4 * mesh.forces().col(_v2);
    F.col(_v9) =
        bu5 * mesh.forces().col(_v1) + bv5 * mesh.forces().col(_v4) + bw5 * mesh.forces().col(_v3);

    F.col(_v5p) = F.col(_v5);
    F.col(_v6p) = F.col(_v6);
    F.col(_v7p) = F.col(_v7);

    /**
     * Transfer velocities to subdivided mesh elements
     */
    V.col(_v5) = (1 - t1) * mesh.velocities().col(_v1) + t1 * mesh.velocities().col(_v4);
    V.col(_v6) = (1 - t2) * mesh.velocities().col(_v2) + t2 * mesh.velocities().col(_v4);
    V.col(_v7) = (1 - t3) * mesh.velocities().col(_v2) + t3 * mesh.velocities().col(_v3);
    V.col(_v8) = bu4 * mesh.velocities().col(_v1) + bv4 * mesh.velocities().col(_v3) +
                 bw4 * mesh.velocities().col(_v2);
    V.col(_v9) = bu5 * mesh.velocities().col(_v1) + bv5 * mesh.velocities().col(_v4) +
                 bw5 * mesh.velocities().col(_v3);

    V.col(_v5p) = V.col(_v5);
    V.col(_v6p) = V.col(_v6);
    V.col(_v7p) = V.col(_v7);

    /**
     * Transfer masses to subdivided mesh elements
     */
    M(_v5) = (1 - t1) * mesh.masses()(_v1) + t1 * mesh.masses()(_v4);
    M(_v6) = (1 - t2) * mesh.masses()(_v2) + t2 * mesh.masses()(_v4);
    M(_v7) = (1 - t3) * mesh.masses()(_v2) + t3 * mesh.masses()(_v3);
    M(_v8) = bu4 * mesh.masses()(_v1) + bv4 * mesh.masses()(_v3) + bw4 * mesh.masses()(_v2);
    M(_v9) = bu5 * mesh.masses()(_v1) + bv5 * mesh.masses()(_v4) + bw5 * mesh.masses()(_v3);

    M(_v5p) = M(_v5);
    M(_v6p) = M(_v6);
    M(_v7p) = M(_v7);

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

    return std::make_tuple(P, T, M, V, F);
}

} // namespace cutting
} // namespace physics
} // namespace sbs