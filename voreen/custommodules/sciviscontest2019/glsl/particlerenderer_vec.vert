#version 330



layout(location = 0) in vec3 in_vertex;
layout(location = 1) in vec3 in_val;

uniform mat4 viewMatrix;

out vec4  geom_position;
out vec3  geom_val;

void main(){
	vec4 pos = viewMatrix* vec4(in_vertex, 1);
	geom_position = pos;
	geom_val = normalize(in_val);
	
	
	gl_Position =  pos;
}