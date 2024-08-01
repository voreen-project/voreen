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

#include "voreen/core/animation/interpolation/stringinterpolationfunctions.h"

namespace voreen {

StringStartInterpolationFunction::StringStartInterpolationFunction() {}

std::string StringStartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string StringStartInterpolationFunction::getCategory() const {
    return "boolean";
}

std::string StringStartInterpolationFunction::interpolate(std::string startvalue, std::string endvalue, float time) const {
    if (time < 1.f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<std::string>* StringStartInterpolationFunction::create() const {
    return new StringStartInterpolationFunction();
}

StringEndInterpolationFunction::StringEndInterpolationFunction() {}

std::string StringEndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}

std::string StringEndInterpolationFunction::getCategory() const {
    return "boolean";
}

std::string StringEndInterpolationFunction::interpolate(std::string startvalue, std::string endvalue, float time) const {
    if (time > 0.f)
        return endvalue;
    else
        return startvalue;
}

InterpolationFunction<std::string>* StringEndInterpolationFunction::create() const {
    return new StringEndInterpolationFunction();
}

StringStartEndInterpolationFunction::StringStartEndInterpolationFunction() {}

std::string StringStartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string StringStartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

std::string StringStartEndInterpolationFunction::interpolate(std::string startvalue, std::string endvalue, float time) const {
    if (time < 0.5f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<std::string>* StringStartEndInterpolationFunction::create() const {
    return new StringStartEndInterpolationFunction();
}

} // namespace voreen
