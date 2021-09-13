#version 330

out uint out_color;

in vec2 frag_coord;
flat in uint frag_instanceID;

void main() {
    if (dot(frag_coord, frag_coord) > 1)
        discard;
    out_color = frag_instanceID;
}
