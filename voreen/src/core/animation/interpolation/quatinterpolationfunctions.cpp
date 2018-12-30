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

#include "voreen/core/animation/interpolation/quatinterpolationfunctions.h"

namespace voreen {

QuatStartInterpolationFunction::QuatStartInterpolationFunction() {}

std::string QuatStartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string QuatStartInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::quat QuatStartInterpolationFunction::interpolate(tgt::quat startvalue, tgt::quat endvalue, float time) const {
    if (time < 1.f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<tgt::quat>* QuatStartInterpolationFunction::create() const {
    return new QuatStartInterpolationFunction();
}

QuatEndInterpolationFunction::QuatEndInterpolationFunction() {}

std::string QuatEndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}

std::string QuatEndInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::quat QuatEndInterpolationFunction::interpolate(tgt::quat startvalue, tgt::quat endvalue, float time) const {
    if (time > 0.f)
        return endvalue;
    else
        return startvalue;
}

InterpolationFunction<tgt::quat>* QuatEndInterpolationFunction::create() const {
        return new QuatEndInterpolationFunction();
}

QuatStartEndInterpolationFunction::QuatStartEndInterpolationFunction() {}

std::string QuatStartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string QuatStartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::quat QuatStartEndInterpolationFunction::interpolate(tgt::quat startvalue, tgt::quat endvalue, float time) const {
    if (time < 0.5f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<tgt::quat>* QuatStartEndInterpolationFunction::create() const {
    return new QuatStartEndInterpolationFunction();
}

QuatLinearInterpolationFunction::QuatLinearInterpolationFunction() {}

std::string QuatLinearInterpolationFunction::getGuiName() const {
    return "linear interpolation (Lerp)";
}

std::string QuatLinearInterpolationFunction::getCategory() const {
    return "linear";
}

tgt::quat QuatLinearInterpolationFunction::interpolate(tgt::quat startvalue, tgt::quat endvalue, float time) const {
    return lerpQuat(startvalue, endvalue, time);
}

InterpolationFunction<tgt::quat>* QuatLinearInterpolationFunction::create() const {
    return new QuatLinearInterpolationFunction();
}

QuatSphericalLinearInterpolationFunction::QuatSphericalLinearInterpolationFunction() {}

std::string QuatSphericalLinearInterpolationFunction::getGuiName() const {
    return "spherical linear interpolation (Slerp)";
}

std::string QuatSphericalLinearInterpolationFunction::getCategory() const {
    return "linear";
}

tgt::quat QuatSphericalLinearInterpolationFunction::interpolate(tgt::quat startvalue, tgt::quat endvalue, float time) const {
    return slerpQuat(startvalue, endvalue, time, false);
}

InterpolationFunction<tgt::quat>* QuatSphericalLinearInterpolationFunction::create() const {
    return new QuatSphericalLinearInterpolationFunction();
}

} // namespace voreen
