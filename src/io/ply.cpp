#include "io/ply.h"

#include "io/tokenize.h"
#include "io/endianness.hpp"

#include <algorithm>
#include <set>
#include <sstream>

namespace sbs {
namespace io {

ply_format_t string_to_format(std::string const& s)
{
    ply_format_t format = ply_format_t::ascii;
    if (s == "ascii")
        format = ply_format_t::ascii;
    if (s == "binary_little_endian")
        format = ply_format_t::binary_little_endian;
    if (s == "binary_big_endian")
        format = ply_format_t::binary_big_endian;
    return format;
}

void write_ply(
    std::filesystem::path const& filepath,
    common::geometry_t const& geometry,
    ply_format_t format)
{
    if (!filepath.has_extension() || filepath.extension() != ".ply")
        return;

    std::ofstream ofs{filepath.c_str(), std::ios::binary};

    if (!ofs.is_open())
        return;

    write_ply(ofs, geometry, format);
}

void write_ply(std::ostream& os, common::geometry_t const& geometry, ply_format_t format)
{
    std::string const color_component_type_str  = "uchar";
    std::string const normal_component_type_str = "float";
    std::string const vertex_component_type_str = "float";
    std::string const uv_component_type_str     = "float";
    std::string const indices_size_type_str     = "uchar";
    std::string const indices_element_type_str  = "int";

    std::ostringstream header_stream{};
    header_stream << "ply\n";

    if (format == ply_format_t::ascii)
        header_stream << "format ascii 1.0\n";
    if (format == ply_format_t::binary_little_endian)
        header_stream << "format binary_little_endian 1.0\n";
    if (format == ply_format_t::binary_big_endian)
        header_stream << "format binary_big_endian 1.0\n";

    common::geometry_t::geometry_type_t const geometry_type = geometry.geometry_type;
    bool const is_tet_mesh      = geometry_type == common::geometry_t::geometry_type_t::tetrahedron;
    bool const is_triangle_mesh = geometry_type == common::geometry_t::geometry_type_t::triangle;

    std::uint32_t const vertex_count = static_cast<std::uint32_t>(geometry.positions.size()) / 3u;
    bool const has_normals           = (geometry.normals.size() == geometry.positions.size());
    bool const has_colors            = (geometry.colors.size() == geometry.positions.size());
    bool const has_uvs = (geometry.uvs.size() == static_cast<std::size_t>(vertex_count * 2u));

    std::uint8_t const mesh_element_list_size =
        is_triangle_mesh ? 3u : is_tet_mesh ? 4u : 3u /* default to 3, but this is wrong */;
    std::uint32_t const mesh_element_count =
        static_cast<std::uint32_t>(geometry.indices.size()) / mesh_element_list_size;

    header_stream << "element vertex " << std::to_string(vertex_count) << "\n"
                  << "property " << vertex_component_type_str << " x\n"
                  << "property " << vertex_component_type_str << " y\n"
                  << "property " << vertex_component_type_str << " z\n";

    if (has_normals)
    {
        header_stream << "property " << normal_component_type_str << " nx\n"
                      << "property " << normal_component_type_str << " ny\n"
                      << "property " << normal_component_type_str << " nz\n";
    }
    if (has_colors)
    {
        header_stream << "property " << color_component_type_str << " r\n"
                      << "property " << color_component_type_str << " g\n"
                      << "property " << color_component_type_str << " b\n";
    }
    if (has_uvs)
    {
        header_stream << "property " << uv_component_type_str << " u\n"
                      << "property " << uv_component_type_str << " v\n";
    }

    header_stream << "element " << (is_triangle_mesh ? "face " : "tet ")
                  << std::to_string(mesh_element_count) << "\n";

    header_stream << "property list " << indices_size_type_str << " " << indices_element_type_str
                  << " indices\n";

    header_stream << "end_header\n";

    std::string const header = header_stream.str();
    os << header;

    if (format == ply_format_t::ascii)
    {
        for (std::size_t i = 0u; i < vertex_count; ++i)
        {
            std::uint32_t const vertex_idx = static_cast<std::uint32_t>(i) * 3u;
            std::uint32_t const normal_idx = static_cast<std::uint32_t>(i) * 3u;
            std::uint32_t const color_idx  = static_cast<std::uint32_t>(i) * 3u;
            std::uint32_t const uv_idx     = static_cast<std::uint32_t>(i) * 2u;

            std::ostringstream oss{};
            float const x = geometry.positions[vertex_idx];
            float const y = geometry.positions[vertex_idx + 1];
            float const z = geometry.positions[vertex_idx + 2];

            oss << std::to_string(x) << " " << std::to_string(y) << " " << std::to_string(z);

            if (has_normals)
            {
                oss << " ";
                float const nx = geometry.normals[normal_idx];
                float const ny = geometry.normals[normal_idx + 1];
                float const nz = geometry.normals[normal_idx + 2];
                oss << std::to_string(nx) << " " << std::to_string(ny) << " " << std::to_string(nz);
            }
            if (has_colors)
            {
                oss << " ";
                std::uint8_t const r = geometry.colors[color_idx];
                std::uint8_t const g = geometry.colors[color_idx + 1];
                std::uint8_t const b = geometry.colors[color_idx + 2];
                oss << std::to_string(r) << " " << std::to_string(g) << " " << std::to_string(b);
            }
            if (has_uvs)
            {
                oss << " ";
                float const u = geometry.uvs[uv_idx];
                float const v = geometry.uvs[uv_idx + 1];
                oss << std::to_string(u) << " " << std::to_string(v);
            }

            oss << "\n";
            os << oss.str();
        }
        for (std::size_t i = 0u; i < mesh_element_count; ++i)
        {
            std::ostringstream oss{};
            std::size_t const element_idx = i * static_cast<std::size_t>(mesh_element_list_size);

            oss << std::to_string(mesh_element_list_size);
            for (std::uint8_t j = 0u; j < mesh_element_list_size; ++j)
            {
                oss << " "
                    << std::to_string(geometry.indices[element_idx + static_cast<std::size_t>(j)]);
            }
            oss << "\n";
            os << oss.str();
        }
    }

    auto const write_binary_data = [&](common::geometry_t const& endian_correct_geometry) {
        // write vertices
        {
            for (std::size_t i = 0u; i < vertex_count; ++i)
            {
                std::uint32_t const vertex_idx = static_cast<std::uint32_t>(i) * 3u;
                std::uint32_t const normal_idx = static_cast<std::uint32_t>(i) * 3u;
                std::uint32_t const color_idx  = static_cast<std::uint32_t>(i) * 3u;
                std::uint32_t const uv_idx     = static_cast<std::uint32_t>(i) * 2u;

                float const x = endian_correct_geometry.positions[vertex_idx];
                float const y = endian_correct_geometry.positions[vertex_idx + 1u];
                float const z = endian_correct_geometry.positions[vertex_idx + 2u];

                os.write(
                    reinterpret_cast<const char*>(&x),
                    static_cast<std::streamsize>(sizeof(x)));
                os.write(
                    reinterpret_cast<const char*>(&y),
                    static_cast<std::streamsize>(sizeof(y)));
                os.write(
                    reinterpret_cast<const char*>(&z),
                    static_cast<std::streamsize>(sizeof(z)));

                if (has_normals)
                {
                    float const nx = endian_correct_geometry.normals[normal_idx];
                    float const ny = endian_correct_geometry.normals[normal_idx + 1u];
                    float const nz = endian_correct_geometry.normals[normal_idx + 2u];

                    os.write(
                        reinterpret_cast<const char*>(&nx),
                        static_cast<std::streamsize>(sizeof(nx)));
                    os.write(
                        reinterpret_cast<const char*>(&ny),
                        static_cast<std::streamsize>(sizeof(ny)));
                    os.write(
                        reinterpret_cast<const char*>(&nz),
                        static_cast<std::streamsize>(sizeof(nz)));
                }
                if (has_colors)
                {
                    std::uint8_t const r = endian_correct_geometry.colors[color_idx];
                    std::uint8_t const g = endian_correct_geometry.colors[color_idx + 1u];
                    std::uint8_t const b = endian_correct_geometry.colors[color_idx + 2u];

                    os.write(
                        reinterpret_cast<const char*>(&r),
                        static_cast<std::streamsize>(sizeof(r)));
                    os.write(
                        reinterpret_cast<const char*>(&g),
                        static_cast<std::streamsize>(sizeof(g)));
                    os.write(
                        reinterpret_cast<const char*>(&b),
                        static_cast<std::streamsize>(sizeof(b)));
                }
                if (has_uvs)
                {
                    float const u = endian_correct_geometry.uvs[uv_idx];
                    float const v = endian_correct_geometry.uvs[uv_idx + 1];

                    os.write(
                        reinterpret_cast<const char*>(&u),
                        static_cast<std::streamsize>(sizeof(u)));
                    os.write(
                        reinterpret_cast<const char*>(&v),
                        static_cast<std::streamsize>(sizeof(v)));
                }
            }
        }

        // write mesh elements
        {
            for (std::size_t i = 0u; i < mesh_element_count; ++i)
            {
                std::size_t const element_idx =
                    i * static_cast<std::size_t>(mesh_element_list_size);

                os.write(
                    reinterpret_cast<const char*>(&mesh_element_list_size),
                    static_cast<std::streamsize>(sizeof(mesh_element_list_size)));
                for (int j = 0; j < 3; ++j)
                {
                    std::int32_t const index = endian_correct_geometry.indices[element_idx + j];
                    os.write(
                        reinterpret_cast<const char*>(&index),
                        static_cast<std::streamsize>(sizeof(index)));
                }
            }
        }
    };

    auto const transform_endianness = [](common::geometry_t const& geometry) {
        common::geometry_t endian_correct_geometry{geometry};

        for (float& coordinate : endian_correct_geometry.positions)
            coordinate = reverse_endianness(coordinate);
        for (float& normal_component : endian_correct_geometry.normals)
            normal_component = reverse_endianness(normal_component);
        for (std::uint8_t& pixel_value : endian_correct_geometry.colors)
            pixel_value = reverse_endianness(pixel_value);
        for (float& uv_coordinate : endian_correct_geometry.uvs)
            uv_coordinate = reverse_endianness(uv_coordinate);
        for (int& index : endian_correct_geometry.indices)
            index = reverse_endianness(index);

        return endian_correct_geometry;
    };

    if (format == ply_format_t::binary_little_endian)
    {
        if (!is_machine_little_endian())
        {
            common::geometry_t const endian_correct_geometry = transform_endianness(geometry);
            write_binary_data(endian_correct_geometry);
        }
        else
        {
            write_binary_data(geometry);
        }
    }

    if (format == ply_format_t::binary_big_endian)
    {
        if (!is_machine_big_endian())
        {
            common::geometry_t const endian_correct_geometry = transform_endianness(geometry);
            write_binary_data(endian_correct_geometry);
        }
        else
        {
            write_binary_data(geometry);
        }
    }
}

std::optional<common::geometry_t> read_ply(std::filesystem::path const& path)
{
    if (!path.has_filename())
        return {};

    if (!path.has_extension() || path.extension() != ".ply")
        return {};

    if (!std::filesystem::exists(path))
        return {};

    std::ifstream fs{path.string(), std::ios::binary};

    if (!fs.is_open())
        return {};

    return read_ply(fs);
}

std::optional<common::geometry_t> read_ply(std::istream& is)
{
    ply_header_description_t description;

    auto const is_valid_property_type = [](std::string const& s) -> bool {
        bool const is_char   = (s == "char");
        bool const is_uchar  = (s == "uchar");
        bool const is_short  = (s == "short");
        bool const is_ushort = (s == "ushort");
        bool const is_int    = (s == "int");
        bool const is_uint   = (s == "uint");
        bool const is_float  = (s == "float");
        bool const is_double = (s == "double");

        // clang-format off
        bool const is_valid = 
            is_char || 
            is_uchar || 
            is_short || 
            is_ushort || 
            is_int || 
            is_uint ||
            is_float || 
            is_double;
        // clang-format on

        return is_valid;
    };

    std::string line;
    bool is_ply = false;
    std::optional<ply_element_t> current_element{};

    /**
     * Read ply header
     */
    while (std::getline(is, line))
    {
        auto const tokens = io::tokenize(line);

        if (!is_ply)
        {
            if (tokens.front() != "ply")
                return {};

            is_ply = true;
            continue;
        }

        if (tokens.empty())
            continue;

        if (tokens.front() == "comment")
            continue;

        if (tokens.front() == "format")
        {
            description.format = string_to_format(tokens.at(1));
            // ignore version
            continue;
        }

        if (tokens.front() == "element")
        {
            if (tokens.size() != 3)
                return {};

            if (current_element.has_value())
                description.elements.push_back(current_element.value());

            ply_element_t element{};
            element.name    = tokens.at(1);
            element.count   = std::stoull(tokens.back());
            current_element = element;

            continue;
        }
        if (tokens.front() == "property")
        {
            bool const are_tokens_of_valid_size = tokens.size() == 3 || tokens.size() == 5;
            if (!are_tokens_of_valid_size)
                return {};

            if (!current_element.has_value())
                return {};

            ply_element_t& element             = current_element.value();
            std::string const type_of_property = tokens.at(1);
            bool const is_list                 = tokens.size() == 5 && type_of_property == "list";

            ply_element_property_t property{};
            property.is_list = is_list;
            property.name    = tokens.back();

            if (!is_list)
            {
                bool const should_process_property = is_valid_property_type(type_of_property);
                if (!should_process_property)
                    return {};

                if (tokens.size() != 3u)
                    return {};

                property.type[0] = type_of_property;
            }
            else
            {
                std::string const type_of_list_size    = tokens.at(2);
                std::string const type_of_list_element = tokens.at(3);
                bool const should_process_property = is_valid_property_type(type_of_list_size) &&
                                                     is_valid_property_type(type_of_list_element);
                if (!should_process_property)
                    return {};

                if (tokens.size() != 5u)
                    return {};

                property.type[0] = type_of_list_size;
                property.type[1] = type_of_list_element;
            }

            element.properties.push_back(property);
        }

        if (tokens.front() == "end_header")
        {
            if (current_element.has_value())
                description.elements.push_back(current_element.value());

            break;
        }
    }

    switch (description.format)
    {
        case ply_format_t::ascii: return read_ply_ascii(is, description);
        case ply_format_t::binary_little_endian: return read_ply_binary(is, description);
        case ply_format_t::binary_big_endian: return read_ply_binary(is, description);
        default: return {};
    }
}

template <class T>
T get_property_value_from_stream(
    std::istream& is,
    std::string const& type,
    bool should_reverse_endianness = false)
{
    std::array<std::byte, 8u> memory{};

    if (type == "char")
    {
        is.read(
            reinterpret_cast<char*>(memory.data()),
            static_cast<std::streamsize>(sizeof(std::int8_t)));
        std::int8_t* memory_ptr = reinterpret_cast<std::int8_t*>(memory.data());
        std::int8_t const value =
            should_reverse_endianness ? io::reverse_endianness(*memory_ptr) : *memory_ptr;
        return static_cast<T>(value);
    }
    if (type == "uchar")
    {
        is.read(
            reinterpret_cast<char*>(memory.data()),
            static_cast<std::streamsize>(sizeof(std::uint8_t)));
        std::uint8_t* memory_ptr = reinterpret_cast<std::uint8_t*>(memory.data());
        std::uint8_t const value =
            should_reverse_endianness ? io::reverse_endianness(*memory_ptr) : *memory_ptr;
        return static_cast<T>(value);
    }
    if (type == "short")
    {
        is.read(
            reinterpret_cast<char*>(memory.data()),
            static_cast<std::streamsize>(sizeof(std::int16_t)));
        std::int16_t* memory_ptr = reinterpret_cast<std::int16_t*>(memory.data());
        std::int16_t const value =
            should_reverse_endianness ? io::reverse_endianness(*memory_ptr) : *memory_ptr;
        return static_cast<T>(value);
    }
    if (type == "ushort")
    {
        is.read(
            reinterpret_cast<char*>(memory.data()),
            static_cast<std::streamsize>(sizeof(std::uint16_t)));
        std::uint16_t* memory_ptr = reinterpret_cast<std::uint16_t*>(memory.data());
        std::uint16_t const value =
            should_reverse_endianness ? io::reverse_endianness(*memory_ptr) : *memory_ptr;
        return static_cast<T>(value);
    }
    if (type == "int")
    {
        is.read(
            reinterpret_cast<char*>(memory.data()),
            static_cast<std::streamsize>(sizeof(std::int32_t)));
        std::int32_t* memory_ptr = reinterpret_cast<std::int32_t*>(memory.data());
        std::int32_t const value =
            should_reverse_endianness ? io::reverse_endianness(*memory_ptr) : *memory_ptr;
        return static_cast<T>(value);
    }
    if (type == "uint")
    {
        is.read(
            reinterpret_cast<char*>(memory.data()),
            static_cast<std::streamsize>(sizeof(std::uint32_t)));
        std::uint32_t* memory_ptr = reinterpret_cast<std::uint32_t*>(memory.data());
        std::uint32_t const value =
            should_reverse_endianness ? io::reverse_endianness(*memory_ptr) : *memory_ptr;
        return static_cast<T>(value);
    }
    if (type == "float")
    {
        is.read(
            reinterpret_cast<char*>(memory.data()),
            static_cast<std::streamsize>(sizeof(float)));
        float* memory_ptr = reinterpret_cast<float*>(memory.data());
        float const value =
            should_reverse_endianness ? io::reverse_endianness(*memory_ptr) : *memory_ptr;
        return static_cast<T>(value);
    }
    if (type == "double")
    {
        is.read(
            reinterpret_cast<char*>(memory.data()),
            static_cast<std::streamsize>(sizeof(double)));
        double* memory_ptr = reinterpret_cast<double*>(memory.data());
        double const value =
            should_reverse_endianness ? io::reverse_endianness(*memory_ptr) : *memory_ptr;
        return static_cast<T>(value);
    }

    return T{};
}

template <class T>
std::vector<T> get_property_list_from_stream(
    std::istream& is,
    std::string const& size_type,
    std::string const& list_element_type,
    bool should_reverse_endianness = false)
{
    std::array<std::byte, 8u> memory{};

    std::size_t list_size =
        get_property_value_from_stream<std::size_t>(is, size_type, should_reverse_endianness);

    std::vector<T> value_list{};
    value_list.reserve(list_size);
    for (std::size_t i = 0u; i < list_size; ++i)
    {
        T const value =
            get_property_value_from_stream<T>(is, list_element_type, should_reverse_endianness);
        value_list.push_back(value);
    }
    return value_list;
}

std::optional<common::geometry_t>
read_ply_ascii(std::istream& is, ply_header_description_t const& description)
{
    common::geometry_t geometry;

    std::string line;

    for (auto const& element : description.elements)
    {
        if (element.name == "vertex")
        {
            std::set<std::string> property_names{};
            for (auto const& property : element.properties)
            {
                property_names.insert(property.name);
            }

            if (property_names.count("x") == 1 && property_names.count("y") == 1 &&
                property_names.count("z") == 1)
            {
                geometry.positions.resize(element.count * 3u);
            }
            if (property_names.count("nx") == 1 && property_names.count("ny") == 1 &&
                property_names.count("nz") == 1)
            {
                geometry.normals.resize(element.count * 3u);
            }
            if (property_names.count("r") == 1 && property_names.count("g") == 1 &&
                property_names.count("b") == 1)
            {
                geometry.colors.resize(element.count * 3u);
            }
            if (property_names.count("u") == 1 && property_names.count("v") == 1)
            {
                geometry.uvs.resize(element.count * 2u);
            }
        }
        if (element.name == "face")
        {
            std::set<std::string> property_names{};
            for (auto const& property : element.properties)
            {
                property_names.insert(property.name);
            }
            if (property_names.count("indices") == 1)
            {
                geometry.indices.resize(element.count * 3u);
            }
            geometry.geometry_type = common::geometry_t::geometry_type_t::triangle;
        }
        if (element.name == "tet")
        {
            std::set<std::string> property_names{};
            for (auto const& property : element.properties)
            {
                property_names.insert(property.name);
            }
            if (property_names.count("indices") == 1)
            {
                geometry.indices.resize(element.count * 4u);
            }
            geometry.geometry_type = common::geometry_t::geometry_type_t::tetrahedron;
        }
    }

    for (auto const& element : description.elements)
    {
        for (std::size_t i = 0; i < element.count; ++i)
        {
            std::getline(is, line);
            if (is.bad())
                return {};

            auto const tokens = io::tokenize(line);
            if (element.name == "vertex")
            {
                std::size_t const position_idx = i * 3u;
                std::size_t const normal_idx   = i * 3u;
                std::size_t const color_idx    = i * 3u;
                std::size_t const uv_idx       = i * 2u;

                for (std::size_t j = 0u; j < element.properties.size(); ++j)
                {
                    ply_element_property_t const& property = element.properties[j];

                    if (property.name == "x")
                    {
                        float const value                     = std::stof(tokens.at(j));
                        geometry.positions[position_idx + 0u] = value;
                    }
                    if (property.name == "y")
                    {
                        float const value                     = std::stof(tokens.at(j));
                        geometry.positions[position_idx + 1u] = value;
                    }
                    if (property.name == "z")
                    {
                        float const value                     = std::stof(tokens.at(j));
                        geometry.positions[position_idx + 2u] = value;
                    }
                    if (property.name == "nx")
                    {
                        float const value                 = std::stof(tokens.at(j));
                        geometry.normals[normal_idx + 0u] = value;
                    }
                    if (property.name == "ny")
                    {
                        float const value                 = std::stof(tokens.at(j));
                        geometry.normals[normal_idx + 1u] = value;
                    }
                    if (property.name == "nz")
                    {
                        float const value                 = std::stof(tokens.at(j));
                        geometry.normals[normal_idx + 2u] = value;
                    }
                    if (property.name == "r")
                    {
                        int const value                 = std::stoi(tokens.at(j));
                        geometry.colors[color_idx + 0u] = static_cast<std::uint8_t>(value);
                    }
                    if (property.name == "g")
                    {
                        int const value                 = std::stoi(tokens.at(j));
                        geometry.colors[color_idx + 1u] = static_cast<std::uint8_t>(value);
                    }
                    if (property.name == "b")
                    {
                        int const value                 = std::stoi(tokens.at(j));
                        geometry.colors[color_idx + 2u] = static_cast<std::uint8_t>(value);
                    }
                    if (property.name == "u")
                    {
                        float const value         = std::stof(tokens.at(j));
                        geometry.uvs[uv_idx + 0u] = value;
                    }
                    if (property.name == "v")
                    {
                        float const value         = std::stof(tokens.at(j));
                        geometry.uvs[uv_idx + 1u] = value;
                    }
                }
            }
            if (element.name == "face")
            {
                std::size_t const element_idx = i * 3u;

                for (std::size_t j = 0u; j < element.properties.size(); ++j)
                {
                    ply_element_property_t const& property = element.properties[j];
                    if (property.name != "indices")
                        continue;

                    if (!property.is_list)
                        continue;

                    int const num_indices_per_element = std::stoi(tokens.at(0));
                    if (num_indices_per_element != 3)
                        continue;

                    int const v1 = std::stoi(tokens.at(1));
                    int const v2 = std::stoi(tokens.at(2));
                    int const v3 = std::stoi(tokens.at(3));

                    geometry.indices[element_idx + 0u] = v1;
                    geometry.indices[element_idx + 1u] = v2;
                    geometry.indices[element_idx + 2u] = v3;
                }
            }
            if (element.name == "tet")
            {
                std::size_t const element_idx = i * 4u;

                for (std::size_t j = 0u; j < element.properties.size(); ++j)
                {
                    ply_element_property_t const& property = element.properties[j];
                    if (property.name != "indices")
                        continue;

                    if (!property.is_list)
                        continue;

                    int const num_indices_per_element = std::stoi(tokens.at(0));
                    if (num_indices_per_element != 4)
                        continue;

                    int const v1 = std::stoi(tokens.at(1));
                    int const v2 = std::stoi(tokens.at(2));
                    int const v3 = std::stoi(tokens.at(3));
                    int const v4 = std::stoi(tokens.at(4));

                    geometry.indices[element_idx + 0u] = v1;
                    geometry.indices[element_idx + 1u] = v2;
                    geometry.indices[element_idx + 2u] = v3;
                    geometry.indices[element_idx + 3u] = v4;
                }
            }
        }
    }

    return geometry;
}

std::optional<common::geometry_t>
read_ply_binary(std::istream& is, ply_header_description_t const& description)
{
    common::geometry_t geometry{};

    for (auto const& element : description.elements)
    {
        if (element.name == "vertex")
        {
            std::size_t property_storage_size = 0u;
            std::set<std::string> property_names{};
            for (auto const& property : element.properties)
            {
                property_names.insert(property.name);
            }

            if (property_names.count("x") == 1 && property_names.count("y") == 1 &&
                property_names.count("z") == 1)
            {
                geometry.positions.resize(element.count * 3u);
            }
            if (property_names.count("nx") == 1 && property_names.count("ny") == 1 &&
                property_names.count("nz") == 1)
            {
                geometry.normals.resize(element.count * 3u);
            }
            if (property_names.count("r") == 1 && property_names.count("g") == 1 &&
                property_names.count("b") == 1)
            {
                geometry.colors.resize(element.count * 3u);
            }
            if (property_names.count("u") == 1 && property_names.count("v") == 1)
            {
                geometry.uvs.resize(element.count * 2u);
            }
        }
        if (element.name == "face")
        {
            std::set<std::string> property_names{};
            for (auto const& property : element.properties)
            {
                property_names.insert(property.name);
            }

            if (property_names.count("indices") == 1)
            {
                geometry.indices.resize(element.count * 3u);
            }
            geometry.geometry_type = common::geometry_t::geometry_type_t::triangle;
        }
        if (element.name == "tet")
        {
            std::set<std::string> property_names{};
            for (auto const& property : element.properties)
            {
                property_names.insert(property.name);
            }

            if (property_names.count("indices") == 1)
            {
                geometry.indices.resize(element.count * 4u);
            }
            geometry.geometry_type = common::geometry_t::geometry_type_t::tetrahedron;
        }
    }

    bool const is_file_little_endian     = description.format == ply_format_t::binary_little_endian;
    bool const is_file_big_endian        = description.format == ply_format_t::binary_big_endian;
    bool const should_reverse_endianness = (io::is_machine_big_endian() && is_file_little_endian) ||
                                           (io::is_machine_little_endian() && is_file_big_endian);

    for (auto const& element : description.elements)
    {
        for (std::size_t i = 0; i < element.count; ++i)
        {
            if (is.bad())
                return {};

            if (element.name == "vertex")
            {
                std::size_t const position_idx = i * 3u;
                std::size_t const normal_idx   = i * 3u;
                std::size_t const color_idx    = i * 3u;
                std::size_t const uv_idx       = i * 2u;

                for (std::size_t j = 0u; j < element.properties.size(); ++j)
                {
                    ply_element_property_t const& property = element.properties[j];

                    if (property.name == "x")
                    {
                        float const value = get_property_value_from_stream<float>(
                            is,
                            property.type[0],
                            should_reverse_endianness);
                        geometry.positions[position_idx + 0u] = value;
                    }
                    if (property.name == "y")
                    {
                        float const value = get_property_value_from_stream<float>(
                            is,
                            property.type[0],
                            should_reverse_endianness);
                        geometry.positions[position_idx + 1u] = value;
                    }
                    if (property.name == "z")
                    {
                        float const value = get_property_value_from_stream<float>(
                            is,
                            property.type[0],
                            should_reverse_endianness);
                        geometry.positions[position_idx + 2u] = value;
                    }
                    if (property.name == "nx")
                    {
                        float const value = get_property_value_from_stream<float>(
                            is,
                            property.type[0],
                            should_reverse_endianness);
                        geometry.normals[normal_idx + 0u] = value;
                    }
                    if (property.name == "ny")
                    {
                        float const value = get_property_value_from_stream<float>(
                            is,
                            property.type[0],
                            should_reverse_endianness);
                        geometry.normals[normal_idx + 1u] = value;
                    }
                    if (property.name == "nz")
                    {
                        float const value = get_property_value_from_stream<float>(
                            is,
                            property.type[0],
                            should_reverse_endianness);
                        geometry.normals[normal_idx + 2u] = value;
                    }
                    if (property.name == "r")
                    {
                        std::uint8_t const value = get_property_value_from_stream<std::uint8_t>(
                            is,
                            property.type[0],
                            should_reverse_endianness);
                        geometry.colors[color_idx + 0u] = value;
                    }
                    if (property.name == "g")
                    {
                        std::uint8_t const value = get_property_value_from_stream<std::uint8_t>(
                            is,
                            property.type[0],
                            should_reverse_endianness);
                        geometry.colors[color_idx + 1u] = value;
                    }
                    if (property.name == "b")
                    {
                        std::uint8_t const value = get_property_value_from_stream<std::uint8_t>(
                            is,
                            property.type[0],
                            should_reverse_endianness);
                        geometry.colors[color_idx + 2u] = value;
                    }
                    if (property.name == "u")
                    {
                        float const value = get_property_value_from_stream<float>(
                            is,
                            property.type[0],
                            should_reverse_endianness);
                        geometry.uvs[uv_idx + 0u] = value;
                    }
                    if (property.name == "v")
                    {
                        float const value = get_property_value_from_stream<float>(
                            is,
                            property.type[0],
                            should_reverse_endianness);
                        geometry.uvs[uv_idx + 1u] = value;
                    }
                }
            }
            if (element.name == "face")
            {
                std::size_t const element_idx = i * 3u;

                for (std::size_t j = 0u; j < element.properties.size(); ++j)
                {
                    ply_element_property_t const& property = element.properties[j];

                    if (!property.is_list)
                    {
                        get_property_value_from_stream<int>(is, property.type[0]);
                        continue;
                    }

                    if (property.name != "indices")
                    {
                        get_property_list_from_stream<int>(is, property.type[0], property.type[1]);
                        continue;
                    }

                    auto const property_list = get_property_list_from_stream<int>(
                        is,
                        property.type[0],
                        property.type[1],
                        should_reverse_endianness);

                    if (property_list.size() != 3)
                        continue;

                    geometry.indices[element_idx + 0u] = property_list[0];
                    geometry.indices[element_idx + 1u] = property_list[1];
                    geometry.indices[element_idx + 2u] = property_list[2];
                }
            }
            if (element.name == "tet")
            {
                std::size_t const element_idx = i * 3u;

                for (std::size_t j = 0u; j < element.properties.size(); ++j)
                {
                    ply_element_property_t const& property = element.properties[j];

                    if (!property.is_list)
                    {
                        get_property_value_from_stream<int>(is, property.type[0]);
                        continue;
                    }

                    if (property.name != "indices")
                    {
                        get_property_list_from_stream<int>(is, property.type[0], property.type[1]);
                        continue;
                    }

                    auto const property_list = get_property_list_from_stream<int>(
                        is,
                        property.type[0],
                        property.type[1],
                        should_reverse_endianness);

                    if (property_list.size() != 4)
                        continue;

                    geometry.indices[element_idx + 0u] = property_list[0];
                    geometry.indices[element_idx + 1u] = property_list[1];
                    geometry.indices[element_idx + 2u] = property_list[2];
                    geometry.indices[element_idx + 3u] = property_list[3];
                }
            }
        }
    }

    if (is.bad())
        return {};

    return geometry;
}

} // namespace io
} // namespace sbs