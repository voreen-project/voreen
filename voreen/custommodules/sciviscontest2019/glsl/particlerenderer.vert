#version 330



layout(location = 0) in vec3 in_vertex;
layout(location = 1) in float in_val;

uniform mat4 viewMatrix;

out vec4  geom_position;
out float  geom_val;

void main(){
	vec4 pos = viewMatrix* vec4(in_vertex, 1);
	geom_position = pos;
	geom_val = in_val;
	
	
	gl_Position =  pos;
}