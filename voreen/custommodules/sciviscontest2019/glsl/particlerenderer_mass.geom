#version 330

layout(points) in;
layout(triangle_strip, max_vertices=4) out;

in vec4  geom_position[1];

out vec2 frag_coord;



uniform mat4 projectionMatrix;
uniform float radius_;


void main() {
	float r = radius_;
	vec4 up = vec4(0, 1, 0, 0);
	vec4 right = vec4(1, 0, 0, 0);

    frag_coord  = vec2(1.0, 1.0);
    gl_Position = projectionMatrix*(r*right+r*up+geom_position[0]);
    gl_Position = gl_Position;
    EmitVertex();

    frag_coord  = vec2(-1.0, 1.0);
    gl_Position = projectionMatrix*(-r*right+r*up+geom_position[0]);
    gl_Position = gl_Position;
    EmitVertex();

    frag_coord  = vec2(1.0, -1.0);
    gl_Position = projectionMatrix*(r*right-r*up+geom_position[0]);
    gl_Position = gl_Position;
    EmitVertex();

    frag_coord  = vec2(-1.0, -1.0);
    gl_Position = projectionMatrix*(-r*right-r*up+geom_position[0]);
    gl_Position = gl_Position;
    EmitVertex(); 

    EndPrimitive();

}
