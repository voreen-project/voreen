#version 330
#define GLSL_VERSION_130
#include "modules/mod_transfunc.frag"

out vec4  out_color;
in  vec2  frag_coord;
in  float frag_bodyValue;


uniform TransFuncParameters transferFunc_;
uniform sampler1D transferFuncTex_;
uniform float alphaFactor_;

uniform float scale_;
uniform float offset_;

void main(){
	if (dot(frag_coord, frag_coord) > 1)
		discard;

	float z= 1- sqrt(frag_coord.x*frag_coord.x + frag_coord.y*frag_coord.y) ;
	vec3 normal = vec3(frag_coord.x,frag_coord.y,z);


	float dcont=max(0.0,dot(normal,vec3(0,0,1)));
	vec3 diffuse=dcont*vec3(0.8,0.8,0.8);
	out_color = vec4(0.5*(applyTF(transferFunc_, transferFuncTex_, frag_bodyValue)).xyz+diffuse,alphaFactor_);
	out_color.a = alphaFactor_;
}
