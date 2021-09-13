#version 330

layout(location = 0) in vec3 in_vertex;
layout(location = 4) in float in_radius;

uniform mat4 viewMatrix;

out float geom_radius;

void main() {
    gl_Position = viewMatrix * vec4(in_vertex, 1);
    geom_radius = in_radius;
}
