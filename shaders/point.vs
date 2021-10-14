#version 330 core

layout (location = 0) in vec3 VertexShaderVertexPosition;
layout (location = 1) in vec3 VertexShaderVertexNormal;
layout (location = 2) in vec3 VertexShaderVertexColor;

out vec3 FragmentShaderVertexPosition;
out vec3 FragmentShaderVertexNormal;
out vec3 FragmentShaderVertexColor;

uniform mat4 view;
uniform mat4 projection;

void main()
{
    FragmentShaderVertexPosition = VertexShaderVertexPosition;
    FragmentShaderVertexNormal = VertexShaderVertexNormal;
    FragmentShaderVertexColor = VertexShaderVertexColor;

    gl_PointSize = 10.0;
    gl_Position = projection * view * vec4(FragmentShaderVertexPosition, 1.0);
}