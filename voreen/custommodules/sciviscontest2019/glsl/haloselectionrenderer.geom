#version 330

layout(points) in;
layout(triangle_strip, max_vertices=4) out;

in vec4  geom_position[1];
in float geom_radius[1];
flat in uint  geom_haloID[1];

out vec2 frag_coord;
flat out uint frag_haloID;


uniform mat4 projectionMatrix;



void main() {
    float r = geom_radius[0];
    vec4 up = vec4(0, 1, 0, 0);
    vec4 right = vec4(1, 0, 0, 0);

    frag_coord  = vec2(1.0, 1.0);
    frag_haloID = geom_haloID[0];
    gl_Position = projectionMatrix*(r*right+r*up+geom_position[0]);
    EmitVertex();

    frag_coord  = vec2(-1.0, 1.0);
    frag_haloID = geom_haloID[0];
    gl_Position = projectionMatrix*(-r*right+r*up+geom_position[0]);
    EmitVertex();

    frag_coord  = vec2(1.0, -1.0);
    frag_haloID = geom_haloID[0];
    gl_Position = projectionMatrix*(r*right-r*up+geom_position[0]);
    EmitVertex();

    frag_coord  = vec2(-1.0, -1.0);
    frag_haloID = geom_haloID[0];
    gl_Position = projectionMatrix*(-r*right-r*up+geom_position[0]);
    EmitVertex();

    EndPrimitive();

}
