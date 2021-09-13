#version 330

layout(points) in;
layout(line_strip, max_vertices=200) out;

in vec4  geom_position[1];
in vec3  geom_vel[1];
in vec3  geom_angular_momenta[1];
in float geom_spin_parameter[1];
in float geom_radius[1];

out float frag_spin;

uniform mat4 projectionMatrix;


void main() {

    frag_spin =geom_spin_parameter[0];
    float r =geom_radius[0];;

    vec3 angularMomenta = geom_angular_momenta[0];	
    vec4 pos = geom_position[0];  //introduce a single vertex at the origin
    vec3 nonCoLinearAngularMomenta = vec3(-angularMomenta.z, angularMomenta.x, angularMomenta.y);
    vec3 sinBase = normalize(cross(angularMomenta, nonCoLinearAngularMomenta))* r*1.2;
    vec3 cosBase = normalize(cross(angularMomenta, sinBase))*r*1.2;

    for(float i = 0; i < 6.38 ; i+=0.1)  //generate vertices at positions on the circumference from 0 to 2*pi
    {
	vec3 newPos= vec3(pos.xyz+cosBase*cos(i)+sinBase*sin(i));
	gl_Position = projectionMatrix* vec4(newPos,1);
	EmitVertex();
    }

    EndPrimitive();


}
