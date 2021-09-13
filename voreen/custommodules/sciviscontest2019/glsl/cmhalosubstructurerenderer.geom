#version 330

layout(points) in;
layout(triangle_strip, max_vertices=4) out;

in vec4  geom_position[1];
in vec3  geom_vel[1];
in float geom_radius[1];
in float geom_mass[1];
flat in int   geom_vertexID[1];

out vec2 frag_coord;
out vec3 frag_vel;
out vec3 frag_view_pos;
out float frag_radius;
out float frag_mass;
flat out int  frag_vertexID;

uniform mat4 projectionMatrix;

void main() {

    float r = geom_radius[0];
    vec4 up = vec4(0, 1, 0, 0);
    vec4 right = vec4(1, 0, 0, 0);

    frag_coord  = vec2(1.0, 1.0);
    frag_vel    = geom_vel[0];
    frag_vertexID = geom_vertexID[0];
    frag_view_pos = geom_position[0].xyz;
    frag_mass = geom_mass[0];
    frag_radius = r;
    gl_Position = projectionMatrix*(r*right+r*up+geom_position[0]);
    EmitVertex();

    frag_coord  = vec2(-1.0, 1.0);
    frag_vel    = geom_vel[0];
    frag_vertexID = geom_vertexID[0];
    frag_view_pos = geom_position[0].xyz;
    frag_mass = geom_mass[0];
    frag_radius = r;
    gl_Position = projectionMatrix*(-r*right+r*up+geom_position[0]);
    EmitVertex();

    frag_coord  = vec2(1.0, -1.0);
    frag_vel    = geom_vel[0];
    frag_vertexID = geom_vertexID[0];
    frag_view_pos = geom_position[0].xyz;
    frag_mass = geom_mass[0];
    frag_radius = r;
    gl_Position = projectionMatrix*(r*right-r*up+geom_position[0]);
    EmitVertex();

    frag_coord  = vec2(-1.0, -1.0);
    frag_vel    = geom_vel[0];
    frag_vertexID = geom_vertexID[0];
    frag_view_pos = geom_position[0].xyz;
    frag_mass = geom_mass[0];
    frag_radius = r;
    gl_Position = projectionMatrix*(-r*right-r*up+geom_position[0]);
    EmitVertex();

    EndPrimitive();

}
