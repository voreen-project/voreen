#version 330
#define GLSL_VERSION_130
#include "modules/mod_transfunc.frag"

out vec4 out_color;

in vec2 frag_coord;


uniform vec4 color_;
uniform float alphaFactor_;


void main(){
	if (dot(frag_coord, frag_coord) > 1)
		discard;
	out_color = color_;
	out_color.a *= alphaFactor_;
}