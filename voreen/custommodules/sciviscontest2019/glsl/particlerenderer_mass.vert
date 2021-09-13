#version 330



layout(location = 0) in vec3 in_vertex;

uniform mat4 viewMatrix;

out vec4  geom_position;

void main(){
	vec4 pos = viewMatrix* vec4(in_vertex, 1);
	geom_position = pos;
	
	
	gl_Position =  pos;
}