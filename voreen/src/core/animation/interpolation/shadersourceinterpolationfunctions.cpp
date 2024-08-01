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

#include "voreen/core/animation/interpolation/shadersourceinterpolationfunctions.h"

namespace voreen {

ShaderSourceStartInterpolationFunction::ShaderSourceStartInterpolationFunction() {}

std::string ShaderSourceStartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string ShaderSourceStartInterpolationFunction::getCategory() const {
    return "boolean";
}

ShaderSource ShaderSourceStartInterpolationFunction::interpolate(ShaderSource startvalue, ShaderSource endvalue, float time) const {
    if (time < 1.f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<ShaderSource>* ShaderSourceStartInterpolationFunction::create() const {
    return new ShaderSourceStartInterpolationFunction();
}

ShaderSourceEndInterpolationFunction::ShaderSourceEndInterpolationFunction() {}

std::string ShaderSourceEndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}

std::string ShaderSourceEndInterpolationFunction::getCategory() const {
    return "boolean";
}

ShaderSource ShaderSourceEndInterpolationFunction::interpolate(ShaderSource startvalue, ShaderSource endvalue, float time) const {
    if (time > 0.f)
        return endvalue;
    else
        return startvalue;
}

InterpolationFunction<ShaderSource>* ShaderSourceEndInterpolationFunction::create() const {
    return new ShaderSourceEndInterpolationFunction();
}

ShaderSourceStartEndInterpolationFunction::ShaderSourceStartEndInterpolationFunction() {}

std::string ShaderSourceStartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string ShaderSourceStartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

ShaderSource ShaderSourceStartEndInterpolationFunction::interpolate(ShaderSource startvalue, ShaderSource endvalue, float time) const {
    if (time < 0.5f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<ShaderSource>* ShaderSourceStartEndInterpolationFunction::create() const {
    return new ShaderSourceStartEndInterpolationFunction();
}

} // namespace voreen
