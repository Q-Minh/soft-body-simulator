#version 330 core

in vec3 FragmentShaderVertexPosition;
in vec3 FragmentShaderVertexNormal;
in vec3 FragmentShaderVertexColor;

out vec4 FragmentColor;

void main()
{
    FragmentColor = vec4(FragmentShaderVertexColor, 1.0f);
}