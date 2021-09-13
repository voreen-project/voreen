#version 330

out vec4 out_color;

in vec2 frag_coord;
in vec3 frag_vel;
in vec3 position;

uniform float alphaFactor_;
uniform vec3 pos;
uniform vec3 light;

uniform float scale_;
uniform float offset_;

void main(){
	if (dot(frag_coord, frag_coord) > 1)
		discard;

        float z= 1- sqrt(frag_coord.x*frag_coord.x + frag_coord.y*frag_coord.y) ;
	vec3 normal = vec3(frag_coord.x,frag_coord.y,z);


        float dcont=max(0.0,dot(normal,vec3(0,0,1)));
	vec3 diffuse=dcont*vec3(0.7,0.7,0.7);
	vec3 color =normalize(frag_vel);
	out_color = vec4(0.5*color+diffuse,alphaFactor_);
}
