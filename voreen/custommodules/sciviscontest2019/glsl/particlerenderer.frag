#version 330
#define GLSL_VERSION_130
#include "modules/mod_transfunc.frag"

out vec4 out_color;

in vec2 frag_coord;
in float frag_val;

uniform TransFuncParameters transferFunc_;
uniform sampler1D transferFuncTex_;
uniform float alphaFactor_;

uniform float scale_;
uniform float offset_;

void main(){
	if (dot(frag_coord, frag_coord) > 1)
		discard;

	out_color = applyTF(transferFunc_, transferFuncTex_, frag_val);
	
	out_color.a *= alphaFactor_;
}