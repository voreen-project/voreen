#version 330

out vec4 out_color;

uniform float alphaFactor_;

void main() {
    out_color = vec4(0.8f,0.8f,0.8f,alphaFactor_);
}
