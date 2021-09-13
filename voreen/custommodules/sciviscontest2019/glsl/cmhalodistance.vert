#version 330


layout(location = 0) in vec3 in_vertex;

uniform mat4 viewMatrix;


void main() {
    gl_Position = viewMatrix * vec4(in_vertex, 1);
}
