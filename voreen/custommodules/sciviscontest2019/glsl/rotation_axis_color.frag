#version 330
#define GLSL_VERSION_130
#include "modules/mod_transfunc.frag"

out vec4 out_color;

in vec2  frag_coord;
in float frag_vel;

uniform TransFuncParameters axisTransferFunc_;
uniform sampler1D axisTransferFuncTex_;
uniform float alphaFactor_;


void main(){
	out_color = applyTF(axisTransferFunc_, axisTransferFuncTex_, frag_vel);
	out_color.a = alphaFactor_ ;
}

