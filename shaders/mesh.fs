#version 330 core

in vec3 FragmentShaderVertexPosition;
in vec3 FragmentShaderVertexNormal;
in vec3 FragmentShaderVertexColor;

out vec4 FragmentColor;

struct DirectionalLightType {
    vec3 direction;
	
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    float exponent;
};

struct PointLightType {
    vec3 position;
    
    float constant;
    float linear;
    float quadratic;
	
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    float exponent;
};

uniform vec3 ViewPosition;
uniform DirectionalLightType DirectionalLight;
uniform PointLightType PointLight;

vec3 GetIntensityFromDirectionalLight(DirectionalLightType Light, vec3 Normal, vec3 ViewingDirection);
vec3 GetIntensityFromPointLight(PointLightType Light, vec3 Normal, vec3 Position, vec3 ViewingDirection);

void main()
{
    vec3 NormalizedVertexNormal = normalize(FragmentShaderVertexNormal);
    vec3 ViewingDirection = normalize(ViewPosition - FragmentShaderVertexPosition);

    vec3 DirectionalLightIntensity = GetIntensityFromDirectionalLight(DirectionalLight, NormalizedVertexNormal, ViewingDirection);
    vec3 PointLightIntensity = GetIntensityFromPointLight(PointLight, NormalizedVertexNormal, FragmentShaderVertexPosition, ViewingDirection);

    vec3 LightIntensity = DirectionalLightIntensity + PointLightIntensity;
    FragmentColor = vec4(LightIntensity, 1.0f);
}

vec3 GetIntensityFromDirectionalLight(DirectionalLightType Light, vec3 Normal, vec3 ViewingDirection)
{
    vec3 LightDirection = normalize(-Light.direction);
    // diffuse shading
    float diff = max(dot(Normal, LightDirection), 0.0);
    // specular shading
    vec3 ReflectDirection = reflect(-LightDirection, Normal);
    float spec = pow(max(dot(ViewingDirection, ReflectDirection), 0.0), Light.exponent);
    // combine results
    vec3 ambient = Light.ambient * FragmentShaderVertexColor;
    vec3 diffuse = Light.diffuse * diff * FragmentShaderVertexColor;
    vec3 specular = Light.specular * spec * FragmentShaderVertexColor;
    return (ambient + diffuse + specular);
}

vec3 GetIntensityFromPointLight(PointLightType Light, vec3 Normal, vec3 Position, vec3 ViewingDirection)
{
    vec3 LightDirection = normalize(Light.position - Position);
    // diffuse shading
    float diff = max(dot(Normal, LightDirection), 0.0);
    // specular shading
    vec3 ReflectDirection = reflect(-LightDirection, Normal);
    float spec = pow(max(dot(ViewingDirection, ReflectDirection), 0.0), Light.exponent);
    // attenuation
    float distance = length(Light.position - Position);
    float attenuation = 1.0 / (Light.constant + Light.linear * distance + Light.quadratic * (distance * distance));    
    // combine results
    vec3 ambient = Light.ambient * FragmentShaderVertexColor;
    vec3 diffuse = Light.diffuse * diff * FragmentShaderVertexColor;
    vec3 specular = Light.specular * spec * FragmentShaderVertexColor;
    ambient *= attenuation;
    diffuse *= attenuation;
    specular *= attenuation;
    return (ambient + diffuse + specular);
}