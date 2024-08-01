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

#include "voreen/core/animation/interpolation/basicfloatinterpolation.h"
//#include <math.h>
#include "tgt/tgt_math.h"

namespace voreen {

float BasicFloatInterpolation::linearInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    return (1-time) * startvalue + time * endvalue;
}

float BasicFloatInterpolation::inQuadInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = pow(time, 2);
    return startvalue + multiplicator * (endvalue-startvalue);
}
float BasicFloatInterpolation::inCubicInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = pow(time, 3);
    return startvalue + multiplicator * (endvalue-startvalue);
}
float BasicFloatInterpolation::inQuartInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = pow(time, 4);
    return startvalue + multiplicator * (endvalue-startvalue);
}
float BasicFloatInterpolation::inQuintInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = pow(time, 5);
    return startvalue + multiplicator * (endvalue-startvalue);
}
float BasicFloatInterpolation::inSineInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = 1.f - (sin((1-time) / (2 * tgt::PIf)));
    return startvalue + multiplicator * (endvalue-startvalue);
}
float BasicFloatInterpolation::inExponentInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = powf(2,10*(time-1)) - 0.0001f;
    return startvalue + multiplicator * (endvalue-startvalue);
}
float BasicFloatInterpolation::inCircInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = 1.f - (sqrt(1.f - pow(time, 2)));
    return startvalue + multiplicator * (endvalue-startvalue);
}

float BasicFloatInterpolation::outQuadInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = 1.f - pow((1-time), 2);
    return startvalue + multiplicator * (endvalue-startvalue);
}

float BasicFloatInterpolation::outCubicInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = 1.f - pow((1-time), 3);
    return startvalue + multiplicator * (endvalue-startvalue);
}

float BasicFloatInterpolation::outQuartInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = 1.f - pow((1-time), 4);
    return startvalue + multiplicator * (endvalue-startvalue);
}

float BasicFloatInterpolation::outQuintInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = 1.f - pow((1-time), 5);
    return startvalue + multiplicator * (endvalue-startvalue);
}

float BasicFloatInterpolation::outSineInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = sin(time / (2.f * tgt::PIf));
    return startvalue + multiplicator * (endvalue-startvalue);
}

float BasicFloatInterpolation::outExponentInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = 1.001f * (1.f - powf(2, -10*time));
    return startvalue + multiplicator * (endvalue-startvalue);
}

float BasicFloatInterpolation::outCircInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    float multiplicator = sqrt(1.f - pow((1-time), 2));
    return startvalue + multiplicator * (endvalue-startvalue);
}

float BasicFloatInterpolation::inOutQuadInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return inQuadInterpolation(startvalue, (startvalue + endvalue) / 2.f, time*2);
    else
        return outQuadInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

float BasicFloatInterpolation::inOutCubicInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return inCubicInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outCubicInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

float BasicFloatInterpolation::inOutQuartInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return inQuartInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outQuartInterpolation((startvalue + endvalue) / 2.f, endvalue , time * 2.f - 1.f);
}

float BasicFloatInterpolation::inOutQuintInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return inQuintInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outQuintInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

float BasicFloatInterpolation::inOutSineInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return inSineInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outSineInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

float BasicFloatInterpolation::inOutExponentInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return inExponentInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outExponentInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

float BasicFloatInterpolation::inOutCircInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return inCircInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outCircInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

float BasicFloatInterpolation::outInQuadInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return outQuadInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inQuadInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

float BasicFloatInterpolation::outInCubicInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return outCubicInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inCubicInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

float BasicFloatInterpolation::outInQuartInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return outQuartInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inQuartInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

float BasicFloatInterpolation::outInQuintInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return outQuintInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inQuintInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

float BasicFloatInterpolation::outInSineInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return outSineInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inSineInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2 - 1.f);
}

float BasicFloatInterpolation::outInExponentInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return outExponentInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inExponentInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

float BasicFloatInterpolation::outInCircInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return startvalue;
    if (time > 1.f)
        return endvalue;
    if (time < 0.5f)
        return outCircInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inCircInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

} // namespace voreen
