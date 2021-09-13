#version 330
#define GLSL_VERSION_130
#include "modules/mod_transfunc.frag"

out vec4 out_color;

in float frag_time;
in vec4 frag_color;

uniform TransFuncParameters transferFunc_;
uniform sampler1D transferFuncTex_;
uniform vec2 interval_;
uniform float alphaFactor_;


void main(){
	if (frag_time < interval_.x  || frag_time > interval_.y) discard;

	float intensity = applyTF(transferFunc_, transferFuncTex_, frag_time).a;
	out_color = frag_color;
	//out_color = applyTF(transferFunc_, transferFuncTex_, frag_time);
	out_color.a *= alphaFactor_ * intensity;
}