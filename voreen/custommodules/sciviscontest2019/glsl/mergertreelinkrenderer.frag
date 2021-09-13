
in vec4 frag_color;
in vec2 frag_normCoords;
in vec2 frag_halfLineDimensions;
in float frag_radiusLeftSquared;
in float frag_radiusRightSquared;
in vec2 frag_texPosSeed;

out vec4 out_color;


uniform sampler2D colorTex_;
uniform float textureZoom_;

float lengthSquared(vec2 v) {
    return dot(v,v);
}
float lengthSquared(vec3 v) {
    return dot(v,v);
}

void main() {
    float normX = abs(frag_normCoords.x);
    float normY = abs(frag_normCoords.y);
    vec2 xshift = vec2(1.0f,0.0f);
    float r = (0.5*(normX*normX)+0.5f);
    if(r<normY
//Cut the halos from the link area. This is necessary for textures with alpha
#ifdef USE_TEXTURES
        || lengthSquared((frag_normCoords+xshift)*frag_halfLineDimensions)<frag_radiusLeftSquared
        || lengthSquared((frag_normCoords-xshift)*frag_halfLineDimensions)<frag_radiusRightSquared
#endif
    ) {
        discard;
    }
#ifdef USE_TEXTURES
    vec2 textureCoords = ((frag_normCoords+frag_texPosSeed)*frag_halfLineDimensions)*textureZoom_;
    vec4 textureColor = texture(colorTex_, textureCoords);
    out_color = textureColor*frag_color;
#else
    float yNormOnRing = normY/r;
    float zNorm = sqrt(1.0f - yNormOnRing*yNormOnRing);
    float zreal = zNorm*r*(frag_halfLineDimensions.y);
    out_color = vec4(frag_color.xyz*zNorm,1);
    gl_FragDepth = 1-zreal*0.1f;
#endif
}
