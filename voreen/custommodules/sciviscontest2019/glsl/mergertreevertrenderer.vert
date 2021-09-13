
#define VERTEX_FLAG_NONE        0u
#define VERTEX_FLAG_SELECTED   (1u<<0)
#define VERTEX_FLAG_MOUSE_OVER (1u<<1)
#define VERTEX_FLAG_ON_PATH    (1u<<2)
#define VERTEX_FLAG_ORPHAN     (1u<<3)

layout(location = 0) in vec3 in_vertex;
layout(location = 1) in uint in_flags;

uniform mat4 viewMatrix;
uniform float radius_;
uniform vec4 colorOrdinary_;
uniform vec4 colorSelected_;
uniform vec4 colorMouseOver_;
uniform vec4 colorOnPath_;
uniform vec4 colorOrphan_;

out vec4 geom_position;
out vec4 geom_color;
out float geom_radius;

void main() {
    vec4 pos = viewMatrix * vec4(in_vertex, 1.0f);
    gl_Position = pos;
    geom_position = pos;
    float r = radius_;
    if((in_flags&VERTEX_FLAG_SELECTED)!=0u) {
        geom_color = colorSelected_;
        r *= 1.5f;
    } else if((in_flags&VERTEX_FLAG_MOUSE_OVER)!=0u) {
        geom_color = colorMouseOver_;
        r *= 1.5f;
    } else if((in_flags&VERTEX_FLAG_ORPHAN)!=0u) {
        geom_color = colorOrphan_;
    } else if((in_flags&VERTEX_FLAG_ON_PATH)!=0u) {
        geom_color = colorOnPath_;
    } else {
        geom_color = colorOrdinary_;
    }
    geom_radius = r;
}
