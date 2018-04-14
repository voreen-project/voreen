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

uniform float randIntensity_;
uniform bool useFilter_;
uniform sampler2D spotNoiseTex_;
uniform vec4 spotColor_;

varying float discard_;

const float PI = 3.14159265358979;
const float SQRT_PI2 = 2.50662827463100;    // == sqrt(2 * PI)
const float SIGMA = 1.41;
const float C = SIGMA / SQRT_PI2;

vec4 spotTexture(vec2 texCoord)     // a circle
{
    texCoord = (texCoord - vec2(0.5)) * 2.0;
    float radius2 = dot(texCoord, texCoord);
    if (radius2 > 1.0)
        return vec4(0.0);

    if (useFilter_ == true)
    {
        float c = 1.0 / C;
        float value = (1.0 - ((c / PI) * exp(-(c * radius2))));
        return vec4(spotColor_.rgb * value, 1.0);
    }
    return spotColor_;
}
/*
vec4 spotTexture(vec2 texCoord)     // a square
{
   return spotColor_;
}
*/
void main() {
    /*if (discard_ > 0.0) {
        discard;
        //gl_FragData[0] = vec4(1.0, 0.0, 0.8, 1.0);
        //gl_FragData[1] = vec4(1.0, 0.0, 0.8, 1.0);
    } else {*/
        vec4 spotColor = spotTexture(gl_TexCoord[2].st);
        vec4 color = texture(spotNoiseTex_, gl_TexCoord[1].st);
        float intensity = randIntensity_;
        if (spotColor.a > 0.0) {
            if (color == vec4(0.0, 0.0, 0.0, 1.0))
                intensity += 0.5;
            color.rgb += (spotColor.rgb * intensity);
        }
        gl_FragData[0] = color;
        gl_FragData[1] = color;
    //}
}
