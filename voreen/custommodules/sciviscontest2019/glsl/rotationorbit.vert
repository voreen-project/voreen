#version 330



layout(location = 0) in vec3 in_vertex;
layout(location = 1) in vec3 in_vel;
layout(location = 2) in vec3 in_angular_momenta;
layout(location = 3) in float in_spin_parameter;
layout(location = 4) in float in_radius;

uniform mat4 viewMatrix;

out vec4  geom_position;
out vec3  geom_vel;
out vec3  geom_angular_momenta;
out float geom_spin_parameter;
out float geom_radius;

void main(){
	vec4 pos = viewMatrix* vec4(in_vertex, 1);
	geom_position = pos;
	geom_vel = in_vel;
	geom_angular_momenta =(viewMatrix* vec4(in_angular_momenta,0)).xyz;
	geom_spin_parameter= in_spin_parameter;
	geom_radius= in_radius;	

	gl_Position =  pos ;
}