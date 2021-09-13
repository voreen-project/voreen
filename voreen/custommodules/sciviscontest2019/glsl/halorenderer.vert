#version 330



layout(location = 0) in vec3 in_vertex;
layout(location = 2) in vec3 in_vel;
layout(location = 4) in float in_radius;
layout(location = 5) in float in_bodyValue;

uniform mat4 viewMatrix;

out vec4  geom_position;
out float geom_bodyValue;
out float geom_radius;
out vec3  geom_vel;


void main(){
	vec4 pos = viewMatrix* vec4(in_vertex, 1);
	geom_position = pos;
	geom_bodyValue = in_bodyValue;
	geom_radius= in_radius;
	geom_vel = in_vel;

	gl_Position =  pos;
}
