/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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
#include "voreen/core/animation/interpolation/mat3interpolationfunctions.h"

namespace voreen {

Mat3StartInterpolationFunction::Mat3StartInterpolationFunction() {}

std::string Mat3StartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string Mat3StartInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::mat3 Mat3StartInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    if (time < 1.f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<tgt::mat3>* Mat3StartInterpolationFunction::create() const {
    return new Mat3StartInterpolationFunction();
}

Mat3EndInterpolationFunction::Mat3EndInterpolationFunction() {}

std::string Mat3EndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}

std::string Mat3EndInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::mat3 Mat3EndInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    if (time > 0.f)
        return endvalue;
    else
        return startvalue;
}

InterpolationFunction<tgt::mat3>* Mat3EndInterpolationFunction::create() const {
    return new Mat3EndInterpolationFunction();
}

Mat3StartEndInterpolationFunction::Mat3StartEndInterpolationFunction() {}

std::string Mat3StartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string Mat3StartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::mat3 Mat3StartEndInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    if (time < 0.5f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<tgt::mat3>* Mat3StartEndInterpolationFunction::create() const {
    return new Mat3StartEndInterpolationFunction();
}

Mat3LinearInterpolationFunction::Mat3LinearInterpolationFunction() {}

std::string Mat3LinearInterpolationFunction::getGuiName() const {
    return "linear interpolation";
}

std::string Mat3LinearInterpolationFunction::getCategory() const {
    return "linear";
}

tgt::mat3 Mat3LinearInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::linearInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3LinearInterpolationFunction::create() const {
    return new Mat3LinearInterpolationFunction();
}

Mat3InQuadInterpolationFunction::Mat3InQuadInterpolationFunction() {}

std::string Mat3InQuadInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat3InQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::mat3 Mat3InQuadInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inQuadInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InQuadInterpolationFunction::create() const {
    return new Mat3InQuadInterpolationFunction();
}

Mat3InCubicInterpolationFunction::Mat3InCubicInterpolationFunction() {}

std::string Mat3InCubicInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat3InCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::mat3 Mat3InCubicInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inCubicInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InCubicInterpolationFunction::create() const {
    return new Mat3InCubicInterpolationFunction();
}

Mat3InQuartInterpolationFunction::Mat3InQuartInterpolationFunction() {}

std::string Mat3InQuartInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat3InQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::mat3 Mat3InQuartInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inQuartInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InQuartInterpolationFunction::create() const {
    return new Mat3InQuartInterpolationFunction();
}

Mat3InQuintInterpolationFunction::Mat3InQuintInterpolationFunction() {}

std::string Mat3InQuintInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat3InQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::mat3 Mat3InQuintInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inQuintInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InQuintInterpolationFunction::create() const {
    return new Mat3InQuintInterpolationFunction();
}

Mat3InSineInterpolationFunction::Mat3InSineInterpolationFunction() {}

std::string Mat3InSineInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat3InSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::mat3 Mat3InSineInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inSineInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InSineInterpolationFunction::create() const {
    return new Mat3InSineInterpolationFunction();
}

Mat3InExponentInterpolationFunction::Mat3InExponentInterpolationFunction() {}

std::string Mat3InExponentInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat3InExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::mat3 Mat3InExponentInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inExponentInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InExponentInterpolationFunction::create() const {
    return new Mat3InExponentInterpolationFunction();
}

Mat3InCircInterpolationFunction::Mat3InCircInterpolationFunction() {}

std::string Mat3InCircInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat3InCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::mat3 Mat3InCircInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inCircInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InCircInterpolationFunction::create() const {
    return new Mat3InCircInterpolationFunction();
}

Mat3OutQuadInterpolationFunction::Mat3OutQuadInterpolationFunction() {}

std::string Mat3OutQuadInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat3OutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::mat3 Mat3OutQuadInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outQuadInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutQuadInterpolationFunction::create() const {
    return new Mat3OutQuadInterpolationFunction();
}

Mat3OutCubicInterpolationFunction::Mat3OutCubicInterpolationFunction() {}

std::string Mat3OutCubicInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat3OutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::mat3 Mat3OutCubicInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outCubicInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutCubicInterpolationFunction::create() const {
    return new Mat3OutCubicInterpolationFunction();
}

Mat3OutQuartInterpolationFunction::Mat3OutQuartInterpolationFunction() {}

std::string Mat3OutQuartInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat3OutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::mat3 Mat3OutQuartInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outQuartInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutQuartInterpolationFunction::create() const {
    return new Mat3OutQuartInterpolationFunction();
}

Mat3OutQuintInterpolationFunction::Mat3OutQuintInterpolationFunction() {}

std::string Mat3OutQuintInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat3OutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::mat3 Mat3OutQuintInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outQuintInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutQuintInterpolationFunction::create() const {
    return new Mat3OutQuintInterpolationFunction();
}

Mat3OutSineInterpolationFunction::Mat3OutSineInterpolationFunction() {}

std::string Mat3OutSineInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat3OutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::mat3 Mat3OutSineInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outSineInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutSineInterpolationFunction::create() const {
    return new Mat3OutSineInterpolationFunction();
}

Mat3OutExponentInterpolationFunction::Mat3OutExponentInterpolationFunction() {}

std::string Mat3OutExponentInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat3OutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::mat3 Mat3OutExponentInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outExponentInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutExponentInterpolationFunction::create() const {
    return new Mat3OutExponentInterpolationFunction();
}

Mat3OutCircInterpolationFunction::Mat3OutCircInterpolationFunction() {}

std::string Mat3OutCircInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat3OutCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::mat3 Mat3OutCircInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outCircInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutCircInterpolationFunction::create() const {
    return new Mat3OutCircInterpolationFunction();
}

Mat3InOutQuadInterpolationFunction::Mat3InOutQuadInterpolationFunction() {}

std::string Mat3InOutQuadInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat3InOutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::mat3 Mat3InOutQuadInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutQuadInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InOutQuadInterpolationFunction::create() const {
    return new Mat3InOutQuadInterpolationFunction();
}

Mat3InOutCubicInterpolationFunction::Mat3InOutCubicInterpolationFunction() {}

std::string Mat3InOutCubicInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat3InOutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::mat3 Mat3InOutCubicInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutCubicInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InOutCubicInterpolationFunction::create() const {
    return new Mat3InOutCubicInterpolationFunction();
}

Mat3InOutQuartInterpolationFunction::Mat3InOutQuartInterpolationFunction() {}

std::string Mat3InOutQuartInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat3InOutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::mat3 Mat3InOutQuartInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutQuartInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InOutQuartInterpolationFunction::create() const {
    return new Mat3InOutQuartInterpolationFunction();
}

Mat3InOutQuintInterpolationFunction::Mat3InOutQuintInterpolationFunction() {}

std::string Mat3InOutQuintInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat3InOutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::mat3 Mat3InOutQuintInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutQuintInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InOutQuintInterpolationFunction::create() const {
    return new Mat3InOutQuintInterpolationFunction();
}

Mat3InOutSineInterpolationFunction::Mat3InOutSineInterpolationFunction() {}

std::string Mat3InOutSineInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat3InOutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::mat3 Mat3InOutSineInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutSineInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InOutSineInterpolationFunction::create() const {
    return new Mat3InOutSineInterpolationFunction();
}

Mat3InOutExponentInterpolationFunction::Mat3InOutExponentInterpolationFunction() {}

std::string Mat3InOutExponentInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat3InOutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::mat3 Mat3InOutExponentInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutExponentInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InOutExponentInterpolationFunction::create() const {
    return new Mat3InOutExponentInterpolationFunction();
}

Mat3InOutCircInterpolationFunction::Mat3InOutCircInterpolationFunction() {}

std::string Mat3InOutCircInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat3InOutCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::mat3 Mat3InOutCircInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutCircInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3InOutCircInterpolationFunction::create() const {
    return new Mat3InOutCircInterpolationFunction();
}

Mat3OutInQuadInterpolationFunction::Mat3OutInQuadInterpolationFunction() {}

std::string Mat3OutInQuadInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat3OutInQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::mat3 Mat3OutInQuadInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInQuadInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutInQuadInterpolationFunction::create() const {
    return new Mat3OutInQuadInterpolationFunction();
}

Mat3OutInCubicInterpolationFunction::Mat3OutInCubicInterpolationFunction() {}

std::string Mat3OutInCubicInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat3OutInCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::mat3 Mat3OutInCubicInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInCubicInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutInCubicInterpolationFunction::create() const {
    return new Mat3OutInCubicInterpolationFunction();
}

Mat3OutInQuartInterpolationFunction::Mat3OutInQuartInterpolationFunction() {}

std::string Mat3OutInQuartInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat3OutInQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::mat3 Mat3OutInQuartInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInQuartInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutInQuartInterpolationFunction::create() const {
    return new Mat3OutInQuartInterpolationFunction();
}

Mat3OutInQuintInterpolationFunction::Mat3OutInQuintInterpolationFunction() {}

std::string Mat3OutInQuintInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat3OutInQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::mat3 Mat3OutInQuintInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInQuintInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutInQuintInterpolationFunction::create() const {
    return new Mat3OutInQuintInterpolationFunction();
}

Mat3OutInSineInterpolationFunction::Mat3OutInSineInterpolationFunction() {}

std::string Mat3OutInSineInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat3OutInSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::mat3 Mat3OutInSineInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInSineInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutInSineInterpolationFunction::create() const {
    return new Mat3OutInSineInterpolationFunction();
}

Mat3OutInExponentInterpolationFunction::Mat3OutInExponentInterpolationFunction() {}

std::string Mat3OutInExponentInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat3OutInExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::mat3 Mat3OutInExponentInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInExponentInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutInExponentInterpolationFunction::create() const {
    return new Mat3OutInExponentInterpolationFunction();
}

Mat3OutInCircInterpolationFunction::Mat3OutInCircInterpolationFunction() {}

std::string Mat3OutInCircInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat3OutInCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::mat3 Mat3OutInCircInterpolationFunction::interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const {
    tgt::mat3 result = tgt::Matrix3<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInCircInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat3>* Mat3OutInCircInterpolationFunction::create() const {
    return new Mat3OutInCircInterpolationFunction();
}

} // namespace voreen
