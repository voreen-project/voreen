#version 330

#define VERTEX_FLAG_NONE        0u
#define VERTEX_FLAG_SELECTED   (1u<<0)
#define VERTEX_FLAG_MOUSE_OVER (1u<<1)
#define VERTEX_FLAG_ON_PATH    (1u<<2)
#define VERTEX_FLAG_ORPHAN     (1u<<3)

#define NO_SELECTION_ID (-1)

layout(location = 0) in vec3 in_vertex;
layout(location = 1) in uint in_flags;
layout(location = 2) in int in_id;

uniform mat4 viewMatrix;
uniform float radius_;

out vec4 geom_position;
out uint geom_instanceID;
out float geom_radius;

void main() {
    vec4 pos = viewMatrix * vec4(in_vertex, 1.0f);
    gl_Position = pos;
    geom_position = pos;
    geom_instanceID = uint(in_id-NO_SELECTION_ID);
    float r = radius_;
    if((in_flags&VERTEX_FLAG_SELECTED)!=0u) {
        r *= 1.5f;
    } else if((in_flags&VERTEX_FLAG_MOUSE_OVER)!=0u) {
        r *= 1.5f;
    }
    geom_radius = r;
}
