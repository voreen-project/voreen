/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

uniform vec2 screenDim_;
uniform vec2 screenDimRCP_;

/*
 * This struct contains information about a texture, like
 * its dimensions. Additionally, the reciprocal
 * values of all parameters are available (suffix RCP) .
 */
struct TextureParameters {
    vec2 dimensions_;        // the texture's resolution, e.g. [256.0, 128.0]
    vec2 dimensionsRCP_;
    mat4 matrix_;            // texture coordinate transformation
    int numChannels_;
    float rwmScale_;         // RealWorldMapping slope
    float rwmOffset_;        // RealWorldMapping intercept
};

// Texture lookup function for 2D textures,
// expecting texture coordinates as pixel coordinates, i.e, [(0,0) , textureSize].
vec4 textureLookup2D(in sampler2D myTexture, in TextureParameters texParams, in vec2 texCoords) {
    vec2 texCoordsNormalized = texCoords * texParams.dimensionsRCP_;
    vec2 texCoordsTransformed = (texParams.matrix_ * vec4(texCoordsNormalized, 0.0, 1.0)).xy;
    return texture(myTexture, texCoordsTransformed);

}

// Texture lookup function for 2D textures,
// expecting normalized texture coordinates, i.e., [0,1].
vec4 textureLookup2Dnormalized(in sampler2D myTexture, in TextureParameters texParams, in vec2 texCoords) {
    vec2 texCoordsTransformed = (texParams.matrix_ * vec4(texCoords, 0.0, 1.0)).xy;
    return texture(myTexture, texCoordsTransformed);
}

// Standard texture lookup function for RenderPort images.
// Texture coordinates are expected in fragment coordinates, i.e, [(0,0) , viewportSize].
//
// @note This function may only be used, if the processor's input and output images (RenderPorts)
//       have the same dimensions.
//
vec4 textureLookup2Dscreen(in sampler2D texture, in TextureParameters texParams, in vec2 texCoords) {
    vec2 texCoordsNormalized = texCoords * screenDimRCP_;
    return textureLookup2Dnormalized(texture, texParams, texCoordsNormalized);
}

// Mimic deprecated texenv function to combine existing fragment color with texture color
vec4 combineColorAndTexture(in vec4 col, in vec4 tex) {
#if defined(TEXTURE_MODE_REPLACE)
    return tex;
#elif defined(TEXTURE_MODE_MODULATE)
    return col * tex;
#elif defined(TEXTURE_MODE_DECAL)
    return vec4((1.0 - tex.a) * col.rgb + tex.a * tex.rgb, col.a);
#elif defined(TEXTURE_MODE_BLEND)
    return vec4((vec3(1.0) - tex.rgb) * col.rgb + col.rgb * tex.rgb, col.a * tex.a);
#else
    // default: return texture color
    return tex;
#endif
}
