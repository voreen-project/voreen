/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#define PI 3.14159265

vec3 fetchLightProbe(in sampler2D lightProbe, in vec3 gradient) {
    // (Dx*r,Dy*r) where r=(1/pi)*acos(Dz)/sqrt(Dx^2 + Dy^2).
    float r = (1.0/PI)*acos(gradient.y)/sqrt(pow(gradient.x, 2.0) + pow(gradient.z, 2.0));
    vec2 texCoord = gradient.xz*r;
    texCoord += 1.0;
    texCoord /= 2.0;
    texCoord.y *= -1.0;
    return texture(lightProbe, texCoord).rgb;
}
