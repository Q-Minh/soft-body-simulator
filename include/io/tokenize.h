#ifndef SBS_IO_TOKENIZE_H
#define SBS_IO_TOKENIZE_H

/**
 * @file
 * @ingroup io
 */

#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace sbs {
namespace io {

/**
 * @ingroup io
 * @brief Tokenize a string by whitespace
 * @param s String to tokenize
 * @return Vector of tokens
 */
inline auto tokenize(std::string const& s) -> std::vector<std::string>
{
    std::istringstream iss(s);
    std::vector<std::string> tokens{std::istream_iterator<std::string>(iss), {}};
    return tokens;
}

} // namespace io
} // namespace sbs

#endif // SBS_IO_TOKENIZE_H