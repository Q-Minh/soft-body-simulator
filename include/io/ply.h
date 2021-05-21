#ifndef SBS_IO_PLY_H
#define SBS_IO_PLY_H

/**
 * @file
 * @ingroup io
 */

#include "io/geometry.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>

namespace sbs {
namespace io {

/**
 * @ingroup io-ply
 * @brief
 * ply format types on disk
 */
enum class ply_format_t { ascii, binary_little_endian, binary_big_endian };

struct ply_element_property_t
{
    std::string name;
    bool is_list = false;
    std::array<std::string, 2> type;
};

struct ply_element_t
{
    std::string name;
    std::size_t count;
    std::vector<ply_element_property_t> properties;
};

struct ply_header_description_t
{
    ply_format_t format = ply_format_t::ascii;
    std::vector<ply_element_t> elements;
};

/**
 * @ingroup io-ply
 * @brief
 * Convert a string representation of a ply format to a ply_format_t enum value
 * @param s The string representation of a ply format
 * @return The corresponding enum value
 */
ply_format_t string_to_format(std::string const& s);

/**
 * @ingroup io-ply
 * @brief
 */
void write_ply(
    std::filesystem::path const& filepath,
    geometry_t const& geometry,
    ply_format_t format = ply_format_t::ascii);

/**
 * @ingroup io-ply
 * @brief
 */
void write_ply(
    std::ostream& os,
    geometry_t const& geometry,
    ply_format_t format = ply_format_t::ascii);

/**
 * @ingroup io-ply
 * @brief
 */
std::optional<geometry_t> read_ply(std::filesystem::path const& path);

/**
 * @ingroup io-ply
 * @brief
 */
std::optional<geometry_t> read_ply(std::istream& is);

/**
 * @ingroup io-ply
 * @brief
 */
std::optional<geometry_t>
read_ply_ascii(std::istream& is, ply_header_description_t const& description);

/**
 * @ingroup io-ply
 * @brief
 */
std::optional<geometry_t>
read_ply_binary(std::istream& is, ply_header_description_t const& description);

} // namespace io
} // namespace sbs

#endif // SBS_IO_PLY_H