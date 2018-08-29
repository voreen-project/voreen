/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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
#include "modules/mod_sampler3d.frag"
#include "modules/mod_transfunc.frag"

// slice/volume
uniform sampler2D sliceTex_;               // slice texture
uniform TextureParameters sliceTexParams_; // slice texture parameters


// transfer functions
uniform TF_SAMPLER_TYPE transFuncTex_;
uniform TransFuncParameters transFuncParams_;

void main() {
    // fetch intensity
    vec4 intensity = textureLookup2Dnormalized(sliceTex_, sliceTexParams_, gl_TexCoord[0].xy);
    intensity *= sliceTexParams_.rwmScale_;
    intensity += sliceTexParams_.rwmOffset_;

    // compositing mode: add channels
    vec4 result;
    result = applyTF(transFuncParams_, transFuncTex_, intensity.a);

    vec4 fragColor = result;
    FragData0 = fragColor;

    if (result.a > 0.0)
        gl_FragDepth = gl_FragCoord.z;
    else
        gl_FragDepth = 1.0;
}

