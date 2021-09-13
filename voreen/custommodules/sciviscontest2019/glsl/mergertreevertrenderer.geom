
layout(points) in;
layout(triangle_strip, max_vertices=4) out;

in vec4  geom_position[1];
in vec4  geom_color[1];
in float geom_radius[1];

out vec4  frag_color;
out vec2  frag_coord;
out vec3  frag_pos;
out float frag_radius;

uniform mat4 projectionMatrix;


void main() {
    float r = geom_radius[0];
    vec4 up = vec4(0, 1, 0, 0);
    vec4 right = vec4(1, 0, 0, 0);

    frag_coord  = vec2(1.0, 1.0);
    frag_color  = geom_color[0];
    gl_Position = projectionMatrix*(r*right+r*up+geom_position[0]);
    frag_pos = gl_Position.xyz;
    frag_radius = r;
    EmitVertex();

    frag_coord  = vec2(-1.0, 1.0);
    frag_color  = geom_color[0];
    gl_Position = projectionMatrix*(-r*right+r*up+geom_position[0]);
    frag_pos = gl_Position.xyz;
    frag_radius = r;
    EmitVertex();

    frag_coord  = vec2(1.0, -1.0);
    frag_color  = geom_color[0];
    gl_Position = projectionMatrix*(r*right-r*up+geom_position[0]);
    frag_pos = gl_Position.xyz;
    frag_radius = r;
    EmitVertex();

    frag_coord  = vec2(-1.0, -1.0);
    frag_color  = geom_color[0];
    gl_Position = projectionMatrix*(-r*right-r*up+geom_position[0]);
    frag_pos = gl_Position.xyz;
    frag_radius = r;
    EmitVertex();

    EndPrimitive();

}
