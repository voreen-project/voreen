/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "modules/mod_transfunc.frag"

in vec2 frag_texcoord;

uniform sampler3D plotData_;
uniform int renderedRuns_[NUM_RUNS]; // NUM_RUNS will be defined by host program.
uniform vec2 rangeX_;
uniform vec2 rangeY_;

// transfer functions
uniform TF_SAMPLER_TYPE transFuncTex_;
uniform TransFuncParameters transFuncParams_;

float mapRange(float valA, float minA, float maxA, float minB, float maxB) {
    return (minB + (maxB - minB) * (valA - minA) / (maxA - minA));
}

void main() {

    float x = mapRange(frag_texcoord.x, 0.0, 1.0, rangeX_.x, rangeX_.y);
    float y = mapRange(frag_texcoord.y, 0.0, 1.0, rangeY_.x, rangeY_.y);

    // Aggregate intensities of all rendered runs.
    float intensity = 0.0;
    int numActive = 0;
    for(int i=0; i < NUM_RUNS; i++) {
        if(renderedRuns_[i] != 0) {
            intensity += texture(plotData_, vec3(x, y, (i+0.5) / NUM_RUNS)).r;
            numActive++;
        }
    }
    
    // Get color from transfer function.
    vec4 fragColor = applyTF(transFuncParams_, transFuncTex_, intensity);
    FragData0 = fragColor;

    if (fragColor.a > 0.0)
        gl_FragDepth = gl_FragCoord.z;
    else
        gl_FragDepth = 1.0;
}

