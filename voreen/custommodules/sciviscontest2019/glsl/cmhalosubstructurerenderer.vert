#version 330



layout(location = 0) in vec3 in_vertex;
layout(location = 1) in vec3 in_vel;
layout(location = 4) in float in_radius;
layout(location = 5) in float in_mass;

uniform mat4 viewMatrix;

out vec4  geom_position;
out vec3  geom_vel;
out float geom_radius;
out float geom_mass;
flat out int   geom_vertexID;


void main(){
    vec4 pos = viewMatrix* vec4(in_vertex, 1);
    geom_position = pos;
    geom_vel = in_vel;
    geom_radius= in_radius;
    geom_vertexID = gl_VertexID;
    geom_mass = in_mass;
    vec3 normVec = in_vertex / abs(in_vertex);
    gl_Position =  pos;
}
