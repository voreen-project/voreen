#version 330

out uint out_haloID;

in vec2 frag_coord;
flat in uint frag_haloID;

void main(){
    if (dot(frag_coord, frag_coord) > 1)
        discard;
    out_haloID = frag_haloID;
}
