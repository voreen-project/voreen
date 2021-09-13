#version 330
layout(lines) in;
layout(line_strip, max_vertices=6) out;

uniform mat4 projectionMatrix;

void main() {

    vec3 first = gl_in[0].gl_Position.xyz;
    vec3 second = gl_in[1].gl_Position.xyz;
    vec3 away = vec3(0,0,1);
    vec3 along = normalize(first - second);
    float l = 0.05*distance(first, second);
    vec3 orth = normalize(cross(away, along))*l;


    gl_Position = projectionMatrix*vec4(first,1);
    EmitVertex();
    gl_Position = projectionMatrix*vec4(second,1);
    EmitVertex();

    EndPrimitive();

    gl_Position = projectionMatrix*vec4(first+orth,1);
    EmitVertex();
    gl_Position = projectionMatrix*vec4(first-orth,1);
    EmitVertex();

    EndPrimitive();

    gl_Position = projectionMatrix*vec4(second+orth,1);
    EmitVertex();
    gl_Position = projectionMatrix*vec4(second-orth,1);
    EmitVertex();

    EndPrimitive();
}
