#ifndef SBS_RENDERING_LIGHT_H
#define SBS_RENDERING_LIGHT_H

namespace sbs {
namespace rendering {

struct ambient_t
{
    ambient_t(float r, float g, float b);
    float r, g, b;
};
struct diffuse_t
{
    diffuse_t(float r, float g, float b);
    float r, g, b;
};
struct specular_t
{
    specular_t(float r, float g, float b, float exp);
    float r, g, b, exp;
};

struct point_light_t
{
    point_light_t(float x, float y, float z);

    float x, y, z;

    struct attenuation_t
    {
        attenuation_t(float constant, float linear, float quadratic);
        float constant, linear, quadratic;
    };

    ambient_t ambient;
    diffuse_t diffuse;
    specular_t specular;
    attenuation_t attenuation;
};

struct directional_light_t
{
    directional_light_t(float dx, float dy, float dz);
    float dx, dy, dz;

    ambient_t ambient;
    diffuse_t diffuse;
    specular_t specular;
};

} // namespace rendering
} // namespace sbs

#endif // SBS_RENDERING_LIGHT_H