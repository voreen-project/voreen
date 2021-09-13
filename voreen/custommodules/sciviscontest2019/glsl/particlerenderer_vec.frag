#version 330

out vec4 out_color;

in vec2 frag_coord;
in vec3 frag_val;

uniform float alphaFactor_;


void main(){
	if (dot(frag_coord, frag_coord) > 1)
		discard;

	out_color = vec4(0.5*frag_val+vec3(0.5), alphaFactor_);
}