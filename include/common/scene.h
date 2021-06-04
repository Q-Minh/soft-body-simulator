#ifndef SBS_COMMON_SCENE_H
#define SBS_COMMON_SCENE_H

#include "node.h"

#include <vector>
#include <memory>

namespace sbs {
namespace common {

struct scene_t
{
    std::vector<std::shared_ptr<node_t>> physics_objects;
    std::vector<std::shared_ptr<node_t>> environment_objects;

    struct ambient_t
    {
        float r, g, b;
    };
    struct diffuse_t
    {
        float r, g, b;
    };
    struct specular_t
    {
        float r, g, b, exp;
    };

    struct point_light_t
    {
        float x, y, z;

        struct attenuation_t
        {
            float constant, linear, quadratic;
        };

        ambient_t ambient;
        diffuse_t diffuse;
        specular_t specular;
        attenuation_t attenuation;
    };

    struct directional_light_t
    {
        float dx, dy, dz;

        ambient_t ambient;
        diffuse_t diffuse;
        specular_t specular;
    };

    directional_light_t directional_light;
    point_light_t point_light;
};

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_SCENE_H