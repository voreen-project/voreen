#version 330

out vec4 out_color;

in vec2 frag_coord;
flat in int frag_val;

uniform float alphaFactor_;


void main(){
	if (dot(frag_coord, frag_coord) > 1)
		discard;
	
	if (frag_val > 3)
		out_color = vec4(vec3(0.5, 0, 0)*(frag_val-4)+vec3(0.5), alphaFactor_);
	else
		out_color = vec4(vec3(0, 0.25, 0)*frag_val+vec3(0.5), alphaFactor_);
}