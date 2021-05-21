#ifndef SBS_IO_ENDIANNESS_HPP
#define SBS_IO_ENDIANNESS_HPP

/**
 * @file
 * @ingroup io
 */

#include <cstddef>
#include <cstdint>

namespace sbs {
namespace io {

/**
 * @ingroup io
 * @brief
 * Checks if current platform is has little endian representation in memory
 * @return True if current plaform is little endian, false otherwise
 */
inline bool is_machine_little_endian()
{
    unsigned int i = 1;
    char* c        = reinterpret_cast<char*>(&i);
    return static_cast<bool>(*c);
}

/**
 * @ingroup io
 * @brief
 * Checks if current platform is has big endian representation in memory
 * @return True if current plaform is big endian, false otherwise
 */
inline bool is_machine_big_endian()
{
    return !is_machine_little_endian();
}

/**
 * @ingroup io
 * @brief
 * Change the endianness of value.
 * @tparam T Type of the value
 * @param value The value for which we wish to change the endianness
 * @return A copy of value with reverse endianness
 */
template <class T>
T reverse_endianness(T value)
{
    union
    {
        T value;
        std::byte u8[sizeof(T)];
    } source, dest;

    source.value = value;
    for (std::size_t i = 0; i < sizeof(T); ++i)
        dest.u8[i] = source.u8[sizeof(T) - i - 1];

    return dest.value;
}

template <>
std::uint8_t reverse_endianness(std::uint8_t value)
{
    return value;
}

template <>
int reverse_endianness(int value)
{
    union
    {
        int value;
        std::byte u8[sizeof(int)];
    } source, dest;
    source.value = value;
    dest.u8[0]   = source.u8[3];
    dest.u8[1]   = source.u8[2];
    dest.u8[2]   = source.u8[1];
    dest.u8[3]   = source.u8[0];

    return dest.value;
}

template <>
float reverse_endianness(float value)
{
    union
    {
        float value;
        std::byte u8[sizeof(float)];
    } source, dest;
    source.value = value;
    dest.u8[0]   = source.u8[3];
    dest.u8[1]   = source.u8[2];
    dest.u8[2]   = source.u8[1];
    dest.u8[3]   = source.u8[0];

    return dest.value;
}

template <>
double reverse_endianness(double value)
{
    union
    {
        double value;
        std::byte u8[sizeof(double)];
    } source, dest;
    source.value = value;
    dest.u8[0]   = source.u8[7];
    dest.u8[1]   = source.u8[6];
    dest.u8[2]   = source.u8[5];
    dest.u8[3]   = source.u8[4];
    dest.u8[4]   = source.u8[3];
    dest.u8[5]   = source.u8[2];
    dest.u8[6]   = source.u8[1];
    dest.u8[7]   = source.u8[0];

    return dest.value;
}

} // namespace io
} // namespace sbs

#endif // SBS_IO_ENDIANNESS_HPP