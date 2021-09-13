
#define VERTEX_FLAG_NONE        0u
#define VERTEX_FLAG_SELECTED   (1u<<0)
#define VERTEX_FLAG_MOUSE_OVER (1u<<1)
#define VERTEX_FLAG_ON_PATH    (1u<<2)
#define VERTEX_FLAG_ORPHAN     (1u<<3)

layout(lines) in;
layout(triangle_strip, max_vertices=4) out;

in vec3 geom_origPos[2];
in uint geom_flags[2];

out vec4 frag_color;
out vec2 frag_normCoords;
out vec2 frag_halfLineDimensions;
out float frag_radiusLeftSquared;
out float frag_radiusRightSquared;
out vec2 frag_texPosSeed;

uniform mat4 projectionMatrix;
uniform float radius_;
uniform float lineWidth_;
uniform vec4 colorOrdinary_;
uniform vec4 colorPath_;

void emitEdgeVertex(vec4 pos, vec2 texPosSeed, vec3 lineSpacingVec, vec2 normcoords, vec2 halfLineDimensions, vec4 color, float rl, float rr) {
    pos.xyz += lineSpacingVec;
    gl_Position = projectionMatrix * pos;
    frag_normCoords = normcoords;
    frag_halfLineDimensions = halfLineDimensions;
    frag_color = color;
    frag_radiusLeftSquared = rl*rl;
    frag_radiusRightSquared = rr*rr;
    frag_texPosSeed = texPosSeed;
    EmitVertex();
}

void main() {
    vec4 geom_color;
    if(((geom_flags[0]&VERTEX_FLAG_ON_PATH)!=0u) && ((geom_flags[1]&VERTEX_FLAG_ON_PATH)!=0u)) {
        geom_color = colorPath_;
    } else {
        geom_color = colorOrdinary_;
    }
    float rl = radius_, rr = radius_;
    if((geom_flags[0]&(VERTEX_FLAG_SELECTED | VERTEX_FLAG_MOUSE_OVER))!=0u) {
        rl *= 1.5f;
    }
    if((geom_flags[1]&(VERTEX_FLAG_SELECTED | VERTEX_FLAG_MOUSE_OVER))!=0u) {
        rr *= 1.5f;
    }
    vec3 away = vec3(0.0f, 0.0f, 1.0f);

    vec3 lineVec = gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz;
    float halfLineLength = length(lineVec.xy)*0.5f;
    float halfLineWidth = lineWidth_*0.5f;
    vec2 halfLineDimensions = vec2(halfLineLength, halfLineWidth);

    vec3 lineWidthDir = normalize(cross(away, lineVec));
    vec3 lineSpacingVec = lineWidthDir*halfLineWidth;
    vec4 pos;

    vec2 texPosSeed = geom_origPos[0].xy;

    emitEdgeVertex(gl_in[0].gl_Position, texPosSeed,  lineSpacingVec, vec2(-1.0f, 1.0f), halfLineDimensions, geom_color, rl, rr);
    emitEdgeVertex(gl_in[0].gl_Position, texPosSeed, -lineSpacingVec, vec2(-1.0f,-1.0f), halfLineDimensions, geom_color, rl, rr);
    emitEdgeVertex(gl_in[1].gl_Position, texPosSeed,  lineSpacingVec, vec2( 1.0f, 1.0f), halfLineDimensions, geom_color, rl, rr);
    emitEdgeVertex(gl_in[1].gl_Position, texPosSeed, -lineSpacingVec, vec2( 1.0f,-1.0f), halfLineDimensions, geom_color, rl, rr);

    EndPrimitive();
}
