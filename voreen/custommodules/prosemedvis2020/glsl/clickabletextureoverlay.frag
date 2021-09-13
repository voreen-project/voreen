
#include "modules/mod_sampler2d.frag"
//include shader libraries (shader modules)
in vec2 frag_texcoord;

uniform sampler2D texture_;

//struct defined in mod_sampler2d.frag containing all needed parameters
uniform TextureParameters textureParameters_;
 
//the color (user-defined by a property)
uniform vec3 color_;

void main() {
    FragData0 = texture(texture_, frag_texcoord)*vec4(color_.r, color_.g, color_.b, 1.0);
}
