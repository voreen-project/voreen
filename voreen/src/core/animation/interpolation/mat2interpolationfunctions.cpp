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
#include "voreen/core/animation/interpolation/mat2interpolationfunctions.h"

namespace voreen {

Mat2StartInterpolationFunction::Mat2StartInterpolationFunction() {}

std::string Mat2StartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string Mat2StartInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::mat2 Mat2StartInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    if (time < 1.f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<tgt::mat2>* Mat2StartInterpolationFunction::create() const {
    return new Mat2StartInterpolationFunction();
}

Mat2EndInterpolationFunction::Mat2EndInterpolationFunction() {}

std::string Mat2EndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}

std::string Mat2EndInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::mat2 Mat2EndInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    if (time > 0.f)
        return endvalue;
    else
        return startvalue;
}

InterpolationFunction<tgt::mat2>* Mat2EndInterpolationFunction::create() const {
    return new Mat2EndInterpolationFunction();
}

Mat2StartEndInterpolationFunction::Mat2StartEndInterpolationFunction() {}

std::string Mat2StartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string Mat2StartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::mat2 Mat2StartEndInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    if (time < 0.5f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<tgt::mat2>* Mat2StartEndInterpolationFunction::create() const {
    return new Mat2StartEndInterpolationFunction();
}

Mat2LinearInterpolationFunction::Mat2LinearInterpolationFunction() {}

std::string Mat2LinearInterpolationFunction::getGuiName() const {
    return "linear interpolation";
}

std::string Mat2LinearInterpolationFunction::getCategory() const {
    return "linear";
}

tgt::mat2 Mat2LinearInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::linearInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2LinearInterpolationFunction::create() const {
    return new Mat2LinearInterpolationFunction();
}

Mat2InQuadInterpolationFunction::Mat2InQuadInterpolationFunction() {}

std::string Mat2InQuadInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat2InQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::mat2 Mat2InQuadInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inQuadInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InQuadInterpolationFunction::create() const {
    return new Mat2InQuadInterpolationFunction();
}

Mat2InCubicInterpolationFunction::Mat2InCubicInterpolationFunction() {}

std::string Mat2InCubicInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat2InCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::mat2 Mat2InCubicInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inCubicInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InCubicInterpolationFunction::create() const {
    return new Mat2InCubicInterpolationFunction();
}

Mat2InQuartInterpolationFunction::Mat2InQuartInterpolationFunction() {}

std::string Mat2InQuartInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat2InQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::mat2 Mat2InQuartInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inQuartInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InQuartInterpolationFunction::create() const {
    return new Mat2InQuartInterpolationFunction();
}

Mat2InQuintInterpolationFunction::Mat2InQuintInterpolationFunction() {}

std::string Mat2InQuintInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat2InQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::mat2 Mat2InQuintInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inQuintInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InQuintInterpolationFunction::create() const {
    return new Mat2InQuintInterpolationFunction();
}

Mat2InSineInterpolationFunction::Mat2InSineInterpolationFunction() {}

std::string Mat2InSineInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat2InSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::mat2 Mat2InSineInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inSineInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InSineInterpolationFunction::create() const {
    return new Mat2InSineInterpolationFunction();
}

Mat2InExponentInterpolationFunction::Mat2InExponentInterpolationFunction() {}

std::string Mat2InExponentInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat2InExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::mat2 Mat2InExponentInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inExponentInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InExponentInterpolationFunction::create() const {
    return new Mat2InExponentInterpolationFunction();
}

Mat2InCircInterpolationFunction::Mat2InCircInterpolationFunction() {}

std::string Mat2InCircInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Mat2InCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::mat2 Mat2InCircInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inCircInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InCircInterpolationFunction::create() const {
    return new Mat2InCircInterpolationFunction();
}

Mat2OutQuadInterpolationFunction::Mat2OutQuadInterpolationFunction() {}

std::string Mat2OutQuadInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat2OutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::mat2 Mat2OutQuadInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outQuadInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutQuadInterpolationFunction::create() const {
    return new Mat2OutQuadInterpolationFunction();
}

Mat2OutCubicInterpolationFunction::Mat2OutCubicInterpolationFunction() {}

std::string Mat2OutCubicInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat2OutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::mat2 Mat2OutCubicInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outCubicInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutCubicInterpolationFunction::create() const {
    return new Mat2OutCubicInterpolationFunction();
}

Mat2OutQuartInterpolationFunction::Mat2OutQuartInterpolationFunction() {}

std::string Mat2OutQuartInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat2OutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::mat2 Mat2OutQuartInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outQuartInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutQuartInterpolationFunction::create() const {
    return new Mat2OutQuartInterpolationFunction();
}

Mat2OutQuintInterpolationFunction::Mat2OutQuintInterpolationFunction() {}

std::string Mat2OutQuintInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat2OutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::mat2 Mat2OutQuintInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outQuintInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutQuintInterpolationFunction::create() const {
    return new Mat2OutQuintInterpolationFunction();
}

Mat2OutSineInterpolationFunction::Mat2OutSineInterpolationFunction() {}

std::string Mat2OutSineInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat2OutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::mat2 Mat2OutSineInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outSineInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutSineInterpolationFunction::create() const {
    return new Mat2OutSineInterpolationFunction();
}

Mat2OutExponentInterpolationFunction::Mat2OutExponentInterpolationFunction() {}

std::string Mat2OutExponentInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat2OutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::mat2 Mat2OutExponentInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outExponentInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutExponentInterpolationFunction::create() const {
    return new Mat2OutExponentInterpolationFunction();
}

Mat2OutCircInterpolationFunction::Mat2OutCircInterpolationFunction() {}

std::string Mat2OutCircInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Mat2OutCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::mat2 Mat2OutCircInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outCircInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutCircInterpolationFunction::create() const {
    return new Mat2OutCircInterpolationFunction();
}

Mat2InOutQuadInterpolationFunction::Mat2InOutQuadInterpolationFunction() {}

std::string Mat2InOutQuadInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat2InOutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::mat2 Mat2InOutQuadInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutQuadInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InOutQuadInterpolationFunction::create() const {
    return new Mat2InOutQuadInterpolationFunction();
}

Mat2InOutCubicInterpolationFunction::Mat2InOutCubicInterpolationFunction() {}

std::string Mat2InOutCubicInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat2InOutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::mat2 Mat2InOutCubicInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutCubicInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InOutCubicInterpolationFunction::create() const {
    return new Mat2InOutCubicInterpolationFunction();
}

Mat2InOutQuartInterpolationFunction::Mat2InOutQuartInterpolationFunction() {}

std::string Mat2InOutQuartInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat2InOutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::mat2 Mat2InOutQuartInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutQuartInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InOutQuartInterpolationFunction::create() const {
    return new Mat2InOutQuartInterpolationFunction();
}

Mat2InOutQuintInterpolationFunction::Mat2InOutQuintInterpolationFunction() {}

std::string Mat2InOutQuintInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat2InOutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::mat2 Mat2InOutQuintInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutQuintInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InOutQuintInterpolationFunction::create() const {
    return new Mat2InOutQuintInterpolationFunction();
}

Mat2InOutSineInterpolationFunction::Mat2InOutSineInterpolationFunction() {}

std::string Mat2InOutSineInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat2InOutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::mat2 Mat2InOutSineInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutSineInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InOutSineInterpolationFunction::create() const {
    return new Mat2InOutSineInterpolationFunction();
}

Mat2InOutExponentInterpolationFunction::Mat2InOutExponentInterpolationFunction() {}

std::string Mat2InOutExponentInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat2InOutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::mat2 Mat2InOutExponentInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutExponentInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InOutExponentInterpolationFunction::create() const {
    return new Mat2InOutExponentInterpolationFunction();
}

Mat2InOutCircInterpolationFunction::Mat2InOutCircInterpolationFunction() {}

std::string Mat2InOutCircInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Mat2InOutCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::mat2 Mat2InOutCircInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::inOutCircInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2InOutCircInterpolationFunction::create() const {
    return new Mat2InOutCircInterpolationFunction();
}

Mat2OutInQuadInterpolationFunction::Mat2OutInQuadInterpolationFunction() {}

std::string Mat2OutInQuadInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat2OutInQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::mat2 Mat2OutInQuadInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInQuadInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutInQuadInterpolationFunction::create() const {
    return new Mat2OutInQuadInterpolationFunction();
}

Mat2OutInCubicInterpolationFunction::Mat2OutInCubicInterpolationFunction() {}

std::string Mat2OutInCubicInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat2OutInCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::mat2 Mat2OutInCubicInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInCubicInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutInCubicInterpolationFunction::create() const {
    return new Mat2OutInCubicInterpolationFunction();
}

Mat2OutInQuartInterpolationFunction::Mat2OutInQuartInterpolationFunction() {}

std::string Mat2OutInQuartInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat2OutInQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::mat2 Mat2OutInQuartInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInQuartInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutInQuartInterpolationFunction::create() const {
    return new Mat2OutInQuartInterpolationFunction();
}

Mat2OutInQuintInterpolationFunction::Mat2OutInQuintInterpolationFunction() {}

std::string Mat2OutInQuintInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat2OutInQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::mat2 Mat2OutInQuintInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInQuintInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutInQuintInterpolationFunction::create() const {
    return new Mat2OutInQuintInterpolationFunction();
}

Mat2OutInSineInterpolationFunction::Mat2OutInSineInterpolationFunction() {}

std::string Mat2OutInSineInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat2OutInSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::mat2 Mat2OutInSineInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInSineInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutInSineInterpolationFunction::create() const {
    return new Mat2OutInSineInterpolationFunction();
}

Mat2OutInExponentInterpolationFunction::Mat2OutInExponentInterpolationFunction() {}

std::string Mat2OutInExponentInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat2OutInExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::mat2 Mat2OutInExponentInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInExponentInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutInExponentInterpolationFunction::create() const {
    return new Mat2OutInExponentInterpolationFunction();
}

Mat2OutInCircInterpolationFunction::Mat2OutInCircInterpolationFunction() {}

std::string Mat2OutInCircInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Mat2OutInCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::mat2 Mat2OutInCircInterpolationFunction::interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const {
    tgt::mat2 result = tgt::Matrix2<float>();
    for (int i = 0; i < result.size; ++i) {
        result.elem[i] = BasicFloatInterpolation::outInCircInterpolation(startvalue.elem[i], endvalue.elem[i], time);
    }
    return result;
}

InterpolationFunction<tgt::mat2>* Mat2OutInCircInterpolationFunction::create() const {
    return new Mat2OutInCircInterpolationFunction();
}

} // namespace voreen
