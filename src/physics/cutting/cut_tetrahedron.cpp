#include "physics/cutting/cut_tetrahedron.h"

namespace sbs {
namespace physics {
namespace cutting {

std::vector<std::pair<Eigen::Matrix3Xd, Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic>>>
cut_tetrahedron(
    Eigen::Matrix3Xd const& V,
    Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic> const& T,
    std::uint32_t tetrahedron,
    std::byte const edge_intersection_mask,
    std::array<Eigen::Vector3d, 6u> const& edge_intersections,
    std::array<Eigen::Vector3d, 4u> const& face_intersections)
{
    auto const v1 = T(0u, tetrahedron);
    auto const v2 = T(1u, tetrahedron);
    auto const v3 = T(2u, tetrahedron);
    auto const v4 = T(3u, tetrahedron);

    auto const& pos1 = V.col(v1);
    auto const& pos2 = V.col(v2);
    auto const& pos3 = V.col(v3);
    auto const& pos4 = V.col(v4);

    tetrahedron_mesh_cutter_t cutter{};
    auto const tets = cutter.subdivide_mesh(
        edge_intersection_mask,
        V,
        T,
        tetrahedron,
        edge_intersections,
        face_intersections);

    return tets;
}

std::vector<tetrahedron_mesh_cutter_t::tetrahedral_mesh_type>
tetrahedron_mesh_cutter_t::subdivide_mesh(
    std::byte const& edge_intersection_mask,
    tetrahedron_mesh_cutter_t::positions_type const& TV,
    tetrahedron_mesh_cutter_t::tetrahedra_type const& TT,
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
        auto const mesh = subdivide_mesh_for_common_case_1(
            TV,
            TT,
            tetrahedron,
            {v1, v3, v4, v2},
            {pe1, pe2, pe3});
        return mesh;
    }
    if (edge_intersection_mask == case_1_134)
    {
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pe2 = edge_intersection_points[e4];
        auto const& pe3 = edge_intersection_points[e3];
        auto const mesh = subdivide_mesh_for_common_case_1(
            TV,
            TT,
            tetrahedron,
            {v2, v4, v3, v1},
            {pe1, pe2, pe3});
        return mesh;
    }
    if (edge_intersection_mask == case_1_236)
    {
        auto const& pe1 = edge_intersection_points[e2];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pe3 = edge_intersection_points[e6];
        auto const mesh = subdivide_mesh_for_common_case_1(
            TV,
            TT,
            tetrahedron,
            {v2, v1, v4, v3},
            {pe1, pe2, pe3});
        return mesh;
    }
    if (edge_intersection_mask == case_1_456)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e5];
        auto const& pe3 = edge_intersection_points[e6];
        auto const mesh = subdivide_mesh_for_common_case_1(
            TV,
            TT,
            tetrahedron,
            {v1, v2, v3, v4},
            {pe1, pe2, pe3});
        return mesh;
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
        auto const mesh = subdivide_mesh_for_common_case_2(
            TV,
            TT,
            tetrahedron,
            {v2, v4, v3, v1},
            {pe1, pe2, pe3, pe4});
        return mesh;
    }
    if (edge_intersection_mask == case_2_1356)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e1];
        auto const& pe3 = edge_intersection_points[e3];
        auto const& pe4 = edge_intersection_points[e6];
        auto const mesh = subdivide_mesh_for_common_case_2(
            TV,
            TT,
            tetrahedron,
            {v2, v3, v1, v4},
            {pe1, pe2, pe3, pe4});
        return mesh;
    }
    if (edge_intersection_mask == case_2_2345)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pe3 = edge_intersection_points[e2];
        auto const& pe4 = edge_intersection_points[e5];
        auto const mesh = subdivide_mesh_for_common_case_2(
            TV,
            TT,
            tetrahedron,
            {v1, v2, v3, v4},
            {pe1, pe2, pe3, pe4});
        return mesh;
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
        auto const mesh = subdivide_mesh_for_common_case_3(
            TV,
            TT,
            tetrahedron,
            {v1, v3, v4, v2},
            {pe1},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_3_2)
    {
        auto const& pe1 = edge_intersection_points[e2];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f4];
        auto const mesh = subdivide_mesh_for_common_case_3(
            TV,
            TT,
            tetrahedron,
            {v3, v4, v1, v2},
            {pe1},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_3_3)
    {
        auto const& pe1 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f4];
        auto const mesh = subdivide_mesh_for_common_case_3(
            TV,
            TT,
            tetrahedron,
            {v1, v4, v2, v3},
            {pe1},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_3_4)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f3];
        auto const mesh = subdivide_mesh_for_common_case_3(
            TV,
            TT,
            tetrahedron,
            {v1, v2, v3, v4},
            {pe1},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_3_5)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f1];
        auto const mesh = subdivide_mesh_for_common_case_3(
            TV,
            TT,
            tetrahedron,
            {v2, v3, v1, v4},
            {pe1},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_3_6)
    {
        auto const& pe1 = edge_intersection_points[e6];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f2];
        auto const mesh = subdivide_mesh_for_common_case_3(
            TV,
            TT,
            tetrahedron,
            {v3, v1, v2, v4},
            {pe1},
            {pf1, pf2});
        return mesh;
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
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v1, v3, v4, v2},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_4_13)
    {
        auto const& pe1 = edge_intersection_points[e3];
        auto const& pe2 = edge_intersection_points[e1];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f3];
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v3, v2, v4, v1},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_4_14)
    {
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pe2 = edge_intersection_points[e4];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f4];
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v2, v4, v3, v1},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_4_15)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e1];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f2];
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v4, v1, v3, v2},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_4_23)
    {
        auto const& pe1 = edge_intersection_points[e2];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f2];
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v2, v1, v4, v3},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_4_25)
    {
        auto const& pe1 = edge_intersection_points[e2];
        auto const& pe2 = edge_intersection_points[e5];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f4];
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v3, v4, v1, v2},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_4_26)
    {
        auto const& pe1 = edge_intersection_points[e6];
        auto const& pe2 = edge_intersection_points[e2];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f3];
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v4, v2, v1, v3},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_4_34)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f1];
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v4, v3, v2, v1},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_4_36)
    {
        auto const& pe1 = edge_intersection_points[e3];
        auto const& pe2 = edge_intersection_points[e6];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f4];
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v1, v4, v2, v3},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_4_45)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e5];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f3];
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v1, v2, v3, v4},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_4_46)
    {
        auto const& pe1 = edge_intersection_points[e6];
        auto const& pe2 = edge_intersection_points[e4];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f2];
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v3, v1, v2, v4},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_4_56)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e6];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f1];
        auto const mesh = subdivide_mesh_for_common_case_4(
            TV,
            TT,
            tetrahedron,
            {v2, v3, v1, v4},
            {pe1, pe2},
            {pf1, pf2});
        return mesh;
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
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v4, v2, v3, v1},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return mesh;
    }
    if (edge_intersection_mask == case_5_126)
    {
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pe2 = edge_intersection_points[e2];
        auto const& pe3 = edge_intersection_points[e6];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f1];
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v1, v3, v4, v2},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_5_135)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e1];
        auto const& pe3 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f3];
        auto const& pf2 = face_intersection_points[f2];
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v4, v1, v3, v2},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_5_136)
    {
        auto const& pe1 = edge_intersection_points[e6];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pe3 = edge_intersection_points[e1];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f2];
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v4, v1, v2, v3},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return mesh;
    }
    if (edge_intersection_mask == case_5_146)
    {
        auto const& pe1 = edge_intersection_points[e1];
        auto const& pe2 = edge_intersection_points[e4];
        auto const& pe3 = edge_intersection_points[e6];
        auto const& pf1 = face_intersection_points[f2];
        auto const& pf2 = face_intersection_points[f4];
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v2, v4, v3, v1},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_5_156)
    {
        auto const& pe1 = edge_intersection_points[e6];
        auto const& pe2 = edge_intersection_points[e5];
        auto const& pe3 = edge_intersection_points[e1];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f3];
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v3, v2, v1, v4},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return mesh;
    }
    if (edge_intersection_mask == case_5_234)
    {
        auto const& pe1 = edge_intersection_points[e2];
        auto const& pe2 = edge_intersection_points[e3];
        auto const& pe3 = edge_intersection_points[e4];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f2];
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v2, v1, v4, v3},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_5_235)
    {
        auto const& pe1 = edge_intersection_points[e3];
        auto const& pe2 = edge_intersection_points[e2];
        auto const& pe3 = edge_intersection_points[e5];
        auto const& pf1 = face_intersection_points[f1];
        auto const& pf2 = face_intersection_points[f3];
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v1, v2, v4, v3},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return mesh;
    }
    if (edge_intersection_mask == case_5_245)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e5];
        auto const& pe3 = edge_intersection_points[e2];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f3];
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v1, v2, v3, v4},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return mesh;
    }
    if (edge_intersection_mask == case_5_246)
    {
        auto const& pe1 = edge_intersection_points[e4];
        auto const& pe2 = edge_intersection_points[e6];
        auto const& pe3 = edge_intersection_points[e2];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f1];
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v1, v3, v2, v4},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return mesh;
    }
    if (edge_intersection_mask == case_5_345)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e4];
        auto const& pe3 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f2];
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v2, v1, v3, v4},
            {pe1, pe2, pe3},
            {pf1, pf2},
            symmetry);
        return mesh;
    }
    if (edge_intersection_mask == case_5_356)
    {
        auto const& pe1 = edge_intersection_points[e5];
        auto const& pe2 = edge_intersection_points[e6];
        auto const& pe3 = edge_intersection_points[e3];
        auto const& pf1 = face_intersection_points[f4];
        auto const& pf2 = face_intersection_points[f1];
        auto const mesh = subdivide_mesh_for_common_case_5(
            TV,
            TT,
            tetrahedron,
            {v2, v3, v1, v4},
            {pe1, pe2, pe3},
            {pf1, pf2});
        return mesh;
    }

    return {};
}

std::vector<tetrahedron_mesh_cutter_t::tetrahedral_mesh_type>
tetrahedron_mesh_cutter_t::subdivide_mesh_for_common_case_1(
    tetrahedron_mesh_cutter_t::positions_type const& TV,
    tetrahedron_mesh_cutter_t::tetrahedra_type const& TT,
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<Eigen::Vector3d, 3u> const& edge_intersection_points)
{
    std::uint32_t const _v1 = TT.col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = TT.col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = TT.col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = TT.col(tetrahedron)(vertex_ordering[3]);

    std::uint32_t constexpr v1 = 0u;
    std::uint32_t constexpr v2 = 1u;
    std::uint32_t constexpr v3 = 2u;
    std::uint32_t constexpr v5 = 3u;
    std::uint32_t constexpr v6 = 4u;
    std::uint32_t constexpr v7 = 5u;

    std::uint32_t constexpr v4  = 0u;
    std::uint32_t constexpr v5p = 1u;
    std::uint32_t constexpr v6p = 2u;
    std::uint32_t constexpr v7p = 3u;

    std::size_t constexpr num_vertices_p1 = 6u;
    std::size_t constexpr num_vertices_p2 = 4u;

    positions_type P1(3u, num_vertices_p1);
    positions_type P2(3u, num_vertices_p2);

    P1.col(v1) = TV.col(_v1);
    P1.col(v2) = TV.col(_v2);
    P1.col(v3) = TV.col(_v3);

    P2.col(v4) = TV.col(_v4);

    P1.col(v5) = edge_intersection_points[0].transpose();
    P1.col(v6) = edge_intersection_points[1].transpose();
    P1.col(v7) = edge_intersection_points[2].transpose();

    P2.col(v5p) = P1.col(v5);
    P2.col(v6p) = P1.col(v6);
    P2.col(v7p) = P1.col(v7);

    std::size_t constexpr num_tetrahedra_t1 = 3u;
    std::size_t constexpr num_tetrahedra_t2 = 1u;
    tetrahedra_type T1(4u, num_tetrahedra_t1);
    tetrahedra_type T2(4u, num_tetrahedra_t2);

    T1.col(0u) = tetrahedron_type{v1, v2, v3, v7};
    T1.col(1u) = tetrahedron_type{v1, v5, v2, v7};
    T1.col(2u) = tetrahedron_type{v6, v5, v7, v2};

    T2.col(0u) = tetrahedron_type{v5p, v6p, v7p, v4};

    std::vector<tetrahedral_mesh_type> tets{};
    tets.reserve(2u);
    tets.push_back({P1, T1});
    tets.push_back({P2, T2});

    return tets;
}

std::vector<tetrahedron_mesh_cutter_t::tetrahedral_mesh_type>
tetrahedron_mesh_cutter_t::subdivide_mesh_for_common_case_2(
    tetrahedron_mesh_cutter_t::positions_type const& TV,
    tetrahedron_mesh_cutter_t::tetrahedra_type const& TT,
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<Eigen::Vector3d, 4u> const& edge_intersection_points)
{
    std::uint32_t const _v1 = TT.col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = TT.col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = TT.col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = TT.col(tetrahedron)(vertex_ordering[3]);

    std::uint32_t constexpr v3 = 0u;
    std::uint32_t constexpr v4 = 1u;
    std::uint32_t constexpr v5 = 2u;
    std::uint32_t constexpr v6 = 3u;
    std::uint32_t constexpr v7 = 4u;
    std::uint32_t constexpr v8 = 5u;

    std::uint32_t constexpr v1  = 0u;
    std::uint32_t constexpr v2  = 1u;
    std::uint32_t constexpr v5p = 2u;
    std::uint32_t constexpr v6p = 3u;
    std::uint32_t constexpr v7p = 4u;
    std::uint32_t constexpr v8p = 5u;

    std::uint32_t constexpr num_vertices_p1 = 6u;
    std::uint32_t constexpr num_vertices_p2 = 6u;
    positions_type P1(3u, num_vertices_p1);
    positions_type P2(3u, num_vertices_p2);

    P2.col(v1) = TV.col(_v1);
    P2.col(v2) = TV.col(_v2);

    P1.col(v3) = TV.col(_v3);
    P1.col(v4) = TV.col(_v4);

    P1.col(v5) = edge_intersection_points[0].transpose();
    P1.col(v6) = edge_intersection_points[1].transpose();
    P1.col(v7) = edge_intersection_points[2].transpose();
    P1.col(v8) = edge_intersection_points[3].transpose();

    P2.col(v5p) = P1.col(v5);
    P2.col(v6p) = P1.col(v6);
    P2.col(v7p) = P1.col(v7);
    P2.col(v8p) = P1.col(v8);

    std::uint32_t constexpr num_tetrahedra_t1 = 3u;
    std::uint32_t constexpr num_tetrahedra_t2 = 3u;
    tetrahedra_type T1(4u, num_tetrahedra_t1);
    tetrahedra_type T2(4u, num_tetrahedra_t2);

    T1.col(0u) = tetrahedron_type{v5, v8, v6, v4};
    T1.col(1u) = tetrahedron_type{v6, v8, v7, v4};
    T1.col(2u) = tetrahedron_type{v6, v7, v3, v4};

    T2.col(0u) = tetrahedron_type{v5p, v6p, v8p, v1};
    T2.col(1u) = tetrahedron_type{v1, v2, v6p, v8p};
    T2.col(2u) = tetrahedron_type{v2, v7p, v6p, v8p};

    std::vector<tetrahedral_mesh_type> tets{};
    tets.reserve(2u);
    tets.push_back({P1, T1});
    tets.push_back({P2, T2});

    return tets;
}

std::vector<tetrahedron_mesh_cutter_t::tetrahedral_mesh_type>
tetrahedron_mesh_cutter_t::subdivide_mesh_for_common_case_3(
    tetrahedron_mesh_cutter_t::positions_type const& TV,
    tetrahedron_mesh_cutter_t::tetrahedra_type const& TT,
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<Eigen::Vector3d, 1u> const& edge_intersection_points,
    std::array<Eigen::Vector3d, 2u> const& face_intersection_points)
{
    std::uint32_t const _v1 = TT.col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = TT.col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = TT.col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = TT.col(tetrahedron)(vertex_ordering[3]);

    std::uint32_t constexpr num_vertices = 8u;
    std::uint32_t constexpr v1           = 0u;
    std::uint32_t constexpr v2           = 1u;
    std::uint32_t constexpr v3           = 2u;
    std::uint32_t constexpr v4           = 3u;
    std::uint32_t constexpr v5           = 4u;
    std::uint32_t constexpr v6           = 5u;
    std::uint32_t constexpr v7           = 6u;
    std::uint32_t constexpr v5p          = 7u;

    positions_type P(3u, num_vertices);

    P.col(v1) = TV.col(_v1);
    P.col(v2) = TV.col(_v2);
    P.col(v3) = TV.col(_v3);
    P.col(v4) = TV.col(_v4);

    P.col(v5) = edge_intersection_points[0].transpose();
    P.col(v6) = face_intersection_points[0].transpose();
    P.col(v7) = face_intersection_points[1].transpose();

    P.col(v5p) = P.col(v5);

    std::uint32_t constexpr num_tetrahedra = 6u;
    tetrahedra_type T(4u, num_tetrahedra);

    T.col(0u) = tetrahedron_type{v5p, v6, v7, v4};
    T.col(1u) = tetrahedron_type{v5, v7, v6, v1};
    T.col(2u) = tetrahedron_type{v6, v3, v7, v4};
    T.col(3u) = tetrahedron_type{v2, v3, v6, v4};
    T.col(4u) = tetrahedron_type{v1, v7, v6, v3};
    T.col(5u) = tetrahedron_type{v1, v2, v3, v6};

    std::vector<tetrahedral_mesh_type> tets{};
    tets.push_back({P, T});
    return tets;
}

std::vector<tetrahedron_mesh_cutter_t::tetrahedral_mesh_type>
tetrahedron_mesh_cutter_t::subdivide_mesh_for_common_case_4(
    tetrahedron_mesh_cutter_t::positions_type const& TV,
    tetrahedron_mesh_cutter_t::tetrahedra_type const& TT,
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<Eigen::Vector3d, 2u> const& edge_intersection_points,
    std::array<Eigen::Vector3d, 2u> const& face_intersection_points)
{
    std::uint32_t const _v1 = TT.col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = TT.col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = TT.col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = TT.col(tetrahedron)(vertex_ordering[3]);

    std::uint32_t constexpr num_vertices = 10u;
    std::uint32_t constexpr v1           = 0u;
    std::uint32_t constexpr v2           = 1u;
    std::uint32_t constexpr v3           = 2u;
    std::uint32_t constexpr v4           = 3u;
    std::uint32_t constexpr v5           = 4u;
    std::uint32_t constexpr v6           = 5u;
    std::uint32_t constexpr v7           = 6u;
    std::uint32_t constexpr v8           = 7u;
    std::uint32_t constexpr v5p          = 8u;
    std::uint32_t constexpr v6p          = 9u;

    positions_type P(3u, num_vertices);

    P.col(v1) = TV.col(_v1);
    P.col(v2) = TV.col(_v2);
    P.col(v3) = TV.col(_v3);
    P.col(v4) = TV.col(_v4);

    P.col(v5) = edge_intersection_points[0].transpose();
    P.col(v6) = edge_intersection_points[1].transpose();
    P.col(v7) = face_intersection_points[0].transpose();
    P.col(v8) = face_intersection_points[1].transpose();

    P.col(v5p) = P.col(v5);
    P.col(v6p) = P.col(v6);

    std::uint32_t constexpr num_tetrahedra = 8u;
    tetrahedra_type T(4u, num_tetrahedra);

    T.col(0u) = tetrahedron_type{v5p, v6p, v7, v4};
    T.col(1u) = tetrahedron_type{v5, v7, v8, v4};
    T.col(2u) = tetrahedron_type{v4, v7, v8, v3};
    T.col(3u) = tetrahedron_type{v7, v6, v5, v2};
    T.col(4u) = tetrahedron_type{v8, v7, v5, v2};
    T.col(5u) = tetrahedron_type{v8, v7, v2, v3};
    T.col(6u) = tetrahedron_type{v5, v8, v2, v1};
    T.col(7u) = tetrahedron_type{v1, v2, v3, v8};

    std::vector<tetrahedral_mesh_type> tets{};
    tets.push_back({P, T});
    return tets;
}

std::vector<tetrahedron_mesh_cutter_t::tetrahedral_mesh_type>
tetrahedron_mesh_cutter_t::subdivide_mesh_for_common_case_5(
    tetrahedron_mesh_cutter_t::positions_type const& TV,
    tetrahedron_mesh_cutter_t::tetrahedra_type const& TT,
    std::uint32_t tetrahedron,
    std::array<std::uint32_t, 4u> const& vertex_ordering,
    std::array<Eigen::Vector3d, 3u> const& edge_intersection_points,
    std::array<Eigen::Vector3d, 2u> const& face_intersection_points,
    bool symmetry)
{
    std::uint32_t const _v1 = TT.col(tetrahedron)(vertex_ordering[0]);
    std::uint32_t const _v2 = TT.col(tetrahedron)(vertex_ordering[1]);
    std::uint32_t const _v3 = TT.col(tetrahedron)(vertex_ordering[2]);
    std::uint32_t const _v4 = TT.col(tetrahedron)(vertex_ordering[3]);

    std::uint32_t constexpr num_vertices = 12u;

    std::uint32_t constexpr v1  = 0u;
    std::uint32_t constexpr v2  = 1u;
    std::uint32_t constexpr v3  = 2u;
    std::uint32_t constexpr v4  = 3u;
    std::uint32_t constexpr v5  = 4u;
    std::uint32_t constexpr v6  = 5u;
    std::uint32_t constexpr v7  = 6u;
    std::uint32_t constexpr v8  = 7u;
    std::uint32_t constexpr v9  = 8u;
    std::uint32_t constexpr v5p = 9u;
    std::uint32_t constexpr v6p = 10u;
    std::uint32_t constexpr v7p = 11u;

    positions_type P(3u, num_vertices);

    P.col(v1) = TV.col(_v1);
    P.col(v2) = TV.col(_v2);
    P.col(v3) = TV.col(_v3);
    P.col(v4) = TV.col(_v4);

    P.col(v5) = edge_intersection_points[0].transpose();
    P.col(v6) = edge_intersection_points[1].transpose();
    P.col(v7) = edge_intersection_points[2].transpose();
    P.col(v8) = face_intersection_points[0].transpose();
    P.col(v9) = face_intersection_points[1].transpose();

    P.col(v5p) = P.col(v5);
    P.col(v6p) = P.col(v6);
    P.col(v7p) = P.col(v7);

    std::uint32_t constexpr num_tetrahedra = 9u;
    tetrahedra_type T(4u, num_tetrahedra);

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

    std::vector<tetrahedral_mesh_type> tets{};
    tets.push_back({P, T});
    return tets;
}

} // namespace cutting
} // namespace physics
} // namespace sbs