#version 330
#define GLSL_VERSION_130
#include "modules/mod_transfunc.frag"

out vec4 out_color;

in float frag_spin;

uniform TransFuncParameters orbitTransferFunc_;
uniform sampler1D orbitTransferFuncTex_;
uniform float alphaFactor_;


void main(){
		out_color = applyTF(orbitTransferFunc_, orbitTransferFuncTex_, frag_spin);
		out_color.a = alphaFactor_;
}

