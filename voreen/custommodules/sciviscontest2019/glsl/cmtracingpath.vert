#version 330



layout(location = 0) in vec3 in_pos;
layout(location = 1) in float in_time;
layout(location = 2) in int in_matrix;
layout(location = 3) in int in_type;


uniform mat4 viewProjectionMatrix_;
uniform mat4 normMats_[101];
uniform vec4 typeColors_[6];

out float frag_time;
out vec4 frag_color;

void main(){
	frag_time = in_time;
	//frag_color = vec4(1.0f,1.0f,1.0f,1.0f);
	frag_color = typeColors_[in_type];
	gl_Position = viewProjectionMatrix_*vec4(in_pos, 1);
}