#include <sbs/rendering/light.h>

namespace sbs {
namespace rendering {

ambient_t::ambient_t(float r, float g, float b) : r(r), g(g), b(b) {}

diffuse_t::diffuse_t(float r, float g, float b) : r(r), g(g), b(b) {}

specular_t::specular_t(float r, float g, float b, float exp) : r(r), g(g), b(b), exp(exp) {}

directional_light_t::directional_light_t(float dx, float dy, float dz)
    : dx(dx),
      dy(dy),
      dz(dz),
      ambient(0.3f, 0.3f, 0.3f),
      diffuse(0.7f, 0.7f, 0.7f),
      specular(0.4f, 0.4f, 0.4f, 100.f)
{
}

point_light_t::attenuation_t::attenuation_t(float constant, float linear, float quadratic)
    : constant(constant), linear(linear), quadratic(quadratic)
{
}

point_light_t::point_light_t(float x, float y, float z)
    : x(x),
      y(y),
      z(z),
      ambient(0.3f, 0.3f, 0.3f),
      diffuse(0.8f, 0.8f, 0.8f),
      specular(0.4f, 0.4f, 0.4f, 100.f),
      attenuation(1.0f, 0.09f, 0.032f)
{
}

} // namespace rendering
} // namespace sbs