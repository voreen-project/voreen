#version 330
layout(lines) in;
layout(triangle_strip, max_vertices=4) out;

in float geom_radius[2];

uniform mat4 projectionMatrix;
uniform float baseWidth;
uniform float spacing;

void emitEdgeVertex(vec3 pos) {
    gl_Position = projectionMatrix * vec4(pos,1);
    EmitVertex();
}

void main() {
    vec3 p0 = gl_in[0].gl_Position.xyz;
    vec3 p1 = gl_in[1].gl_Position.xyz;

    vec3 alongLine = normalize(p1 - p0);
    vec3 orthLine = normalize(cross(alongLine, vec3(0,0,1)))*baseWidth;

    emitEdgeVertex(p0+orthLine+alongLine*(geom_radius[0]*spacing));
    emitEdgeVertex(p0-orthLine+alongLine*(geom_radius[0]*spacing));
    emitEdgeVertex(p1         -alongLine*(geom_radius[1]*spacing));

    EndPrimitive();
}
