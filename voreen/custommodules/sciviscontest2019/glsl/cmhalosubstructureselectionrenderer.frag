#version 330

out uint out_haloID;

in vec2 frag_coord;
flat in uint frag_haloID;
flat in int  frag_vertexID;
in float frag_pos_z;
in float frag_radius;

uniform float hostHaloRingThickness_;
uniform mat4 projectionMatrix;

void main(){
    float oneMinusThickness = 1-hostHaloRingThickness_;
    float fdot = dot(frag_coord, frag_coord);
    if (fdot > 1 || (frag_vertexID==0 && fdot < (oneMinusThickness*oneMinusThickness)))
        discard;
    float dz= sqrt(1 - frag_coord.x*frag_coord.x - frag_coord.y*frag_coord.y);
    vec4 view_z = vec4(0,0,frag_pos_z+dz*frag_radius,1);
    vec4 predevide_z = (projectionMatrix * view_z);

    gl_FragDepth = (0.5*predevide_z.z/predevide_z.w) + 0.5;

    out_haloID = frag_haloID;
}
