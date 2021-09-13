
layout(location = 0) in vec3 in_vertex;
layout(location = 1) in uint in_flags;

uniform mat4 viewMatrix;

out vec3 geom_origPos;
out uint geom_flags;

void main() {
    geom_origPos = in_vertex;
    gl_Position = viewMatrix * vec4(in_vertex, 1);
    geom_flags = in_flags;
}
