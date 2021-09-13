#version 330

#define NO_SELECTION_ID (-1)

layout(location = 0) in vec3 in_vertex;
layout(location = 1) in float in_radius;
layout(location = 2) in int in_haloID;

uniform mat4 viewMatrix;

out vec4  geom_position;
out float geom_radius;
flat out uint  geom_haloID;

void main(){
    vec4 pos = viewMatrix* vec4(in_vertex, 1);
    geom_position = pos;
    geom_radius = in_radius;
    geom_haloID = uint(in_haloID-NO_SELECTION_ID);

    gl_Position =  pos;
}
