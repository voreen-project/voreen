/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
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

#include "modules/mod_sampler2d.frag"

uniform sampler2D colorTex_;
uniform sampler2D depthTex_;
uniform TextureParameters texParams_;

uniform float startThreshold_;

/***
 * The main method.
 ***/
void main() {
    vec4 fragCoord = gl_FragCoord;
    vec2 gradient = textureLookup2Dscreen(colorTex_, texParams_, fragCoord.xy).xy;

    float gradientLength = length(gradient);
    gradientLength = textureLookup2Dscreen(colorTex_, texParams_, fragCoord.xy).b;
    float gradientAngle = atan(gradient.y/gradient.x);
    vec2 offset = vec2(round(cos(gradientAngle)), round(sin(gradientAngle)));

    float myDepth = textureLookup2Dscreen(depthTex_, texParams_, fragCoord.xy).z;
    float forwardNeighbor = textureLookup2Dscreen(colorTex_, texParams_, fragCoord.xy+offset).b;
    float backwardNeighbor = textureLookup2Dscreen(colorTex_, texParams_, fragCoord.xy-offset).b;
    float forwardNeighborDepth = textureLookup2Dscreen(depthTex_, texParams_, fragCoord.xy+offset).z;
    float backwardNeighborDepth = textureLookup2Dscreen(depthTex_, texParams_, fragCoord.xy-offset).z;
    float isMaximum = 0.0;
    if (forwardNeighbor == 0.0) forwardNeighborDepth = -1.0;
    if (backwardNeighbor == 0.0) backwardNeighborDepth = -1.0;
    if (forwardNeighbor < gradientLength && backwardNeighbor < gradientLength) {// non-maximum supression
        //if (forwardNeighborDepth < myDepth && backwardNeighborDepth < myDepth) // non-maximum supression
        if (gradientLength > startThreshold_) // thresholding
            isMaximum = 1.0;
    }

    // now we don't want to go in gradient direction, but along the edge
    offset.x = -offset.x;

    // allow to store also negative values
    offset /= 2.0;
    offset += 0.5;

    if (isMaximum == 1.0)
        FragData0 = vec4(1.0);
    else FragData0 = vec4(0.0);
    gl_FragDepth = myDepth;
}
