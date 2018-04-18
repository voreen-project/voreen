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

#include "voreen/core/animation/interpolation/basicintinterpolation.h"
#include <math.h>
#include "tgt/tgt_math.h"

namespace voreen {

int BasicIntInterpolation::linearInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    return static_cast<int>((1-time) * startvalue + time * endvalue);
}

int BasicIntInterpolation::inQuadInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = pow(time, 2);
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::inCubicInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = pow(time, 3);
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::inQuartInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = pow(time, 4);
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::inQuintInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = pow(time, 5);
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::inSineInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = 1.f - sin((1.f - time) / (2.f * tgt::PIf));
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::inExponentInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = pow(2.f, 10.f * (time - 1.f)) - 0.0001f;
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::inCircInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = 1.f - sqrt(1.f - pow(time, 2));
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::outQuadInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = 1.f - pow((1.f - time), 2);
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::outCubicInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = 1.f - pow((1.f - time), 3);
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::outQuartInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = 1.f - pow((1.f - time), 4);
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::outQuintInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = 1.f - pow((1.f - time), 5);
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::outSineInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = sin(time / (2.f * tgt::PIf));
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::outExponentInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = 1.001f * (1.f - pow(2.f, -10.f * time));
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::outCircInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    float multiplier = sqrt(1.f - pow((1.f - time), 2));
    return static_cast<int>(startvalue + multiplier * (endvalue-startvalue));
}

int BasicIntInterpolation::inOutQuadInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return inQuadInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outQuadInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::inOutCubicInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return inCubicInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outCubicInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::inOutQuartInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return inQuartInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outQuartInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::inOutQuintInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f) return static_cast<int>(startvalue);
    if (time > 1.f) return static_cast<int>(endvalue);
    if (time < 0.5f)
        return inQuintInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outQuintInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::inOutSineInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return inSineInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outSineInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::inOutExponentInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return inExponentInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outExponentInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::inOutCircInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return inCircInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return outCircInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::outInQuadInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return outQuadInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inQuadInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::outInCubicInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return outCubicInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inCubicInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::outInQuartInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return outQuartInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inQuartInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::outInQuintInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return outQuintInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inQuintInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::outInSineInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return outSineInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inSineInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::outInExponentInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return outExponentInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inExponentInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

int BasicIntInterpolation::outInCircInterpolation(float startvalue, float endvalue, float time) {
    if (time < 0.f)
        return static_cast<int>(startvalue);
    if (time > 1.f)
        return static_cast<int>(endvalue);
    if (time < 0.5f)
        return outCircInterpolation(startvalue, (startvalue + endvalue) / 2.f, time * 2.f);
    else
        return inCircInterpolation((startvalue + endvalue) / 2.f, endvalue, time * 2.f - 1.f);
}

} // namespace voreen
