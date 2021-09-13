
out vec4 out_color;

in vec4 frag_color;
in vec2 frag_coord;
in vec3 frag_pos;
in float frag_radius;

uniform sampler2D colorTex_;

void main() {
    if (dot(frag_coord, frag_coord) > 1)
        discard;
#ifdef USE_TEXTURES
    vec4 textureColor = texture(colorTex_, 0.5f*(vec2(1.0f, 1.0f)+frag_coord));
    out_color = textureColor.x*frag_color;
#else
    float dz= sqrt(1 - frag_coord.x*frag_coord.x - frag_coord.y*frag_coord.y);
    float dcont = max(0.0,dz);
    out_color = vec4(pow(dcont,1)*frag_color.xyz, 1);
    gl_FragDepth = 1-(dz*frag_radius*0.1f);
#endif
}
