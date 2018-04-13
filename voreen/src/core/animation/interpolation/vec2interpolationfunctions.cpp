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

#include "voreen/core/animation/interpolation/basicfloatinterpolation.h"
#include "voreen/core/animation/interpolation/vec2interpolationfunctions.h"

namespace voreen {

Vec2StartInterpolationFunction::Vec2StartInterpolationFunction() {}

std::string Vec2StartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string Vec2StartInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::vec2 Vec2StartInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    if (time < 1.f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<tgt::vec2>* Vec2StartInterpolationFunction::create() const {
    return new Vec2StartInterpolationFunction();
}

Vec2EndInterpolationFunction::Vec2EndInterpolationFunction() {}

std::string Vec2EndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}

std::string Vec2EndInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::vec2 Vec2EndInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    if (time > 0.f)
        return endvalue;
    else
        return startvalue;
}

InterpolationFunction<tgt::vec2>* Vec2EndInterpolationFunction::create() const {
    return new Vec2EndInterpolationFunction();
}

Vec2StartEndInterpolationFunction::Vec2StartEndInterpolationFunction() {}

std::string Vec2StartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string Vec2StartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::vec2 Vec2StartEndInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    if (time < 0.5f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<tgt::vec2>* Vec2StartEndInterpolationFunction::create() const {
    return new Vec2StartEndInterpolationFunction();
}

Vec2LinearInterpolationFunction::Vec2LinearInterpolationFunction() {}

std::string Vec2LinearInterpolationFunction::getGuiName() const {
    return "linear interpolation";
}

std::string Vec2LinearInterpolationFunction::getCategory() const {
    return "linear";
}

tgt::vec2 Vec2LinearInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::linearInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::linearInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2LinearInterpolationFunction::create() const {
    return new Vec2LinearInterpolationFunction();
}

Vec2SphericalLinearInterpolationFunction::Vec2SphericalLinearInterpolationFunction() {}

std::string Vec2SphericalLinearInterpolationFunction::getGuiName() const {
        return "spherical linear interpolation";
}

std::string Vec2SphericalLinearInterpolationFunction::getCategory() const {
        return "linear";
}

tgt::vec2 Vec2SphericalLinearInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
        float eps = 0.0001f;
        float dotVal = dot(normalize(startvalue), normalize(endvalue));
        if (dotVal > 0.995f){
            // small angle => linear interpolation
            float compX = BasicFloatInterpolation::linearInterpolation(startvalue.x, endvalue.x, time);
            float compY = BasicFloatInterpolation::linearInterpolation(startvalue.y, endvalue.y, time);
            return tgt::vec2(compX, compY);
        }
        float theta = std::acos(dotVal);
        float compX = (sin((1.f - time) * theta)/(sin(theta) + eps)) * startvalue.x
                      + (sin(time * theta)/(sin(theta) + eps)) * endvalue.x;
        float compY = (sin((1.f - time) * theta)/(sin(theta) + eps)) * startvalue.y
                      + (sin(time * theta)/(sin(theta) + eps)) * endvalue.y;
        return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2SphericalLinearInterpolationFunction::create() const {
        return new Vec2SphericalLinearInterpolationFunction();
}

Vec2InQuadInterpolationFunction::Vec2InQuadInterpolationFunction() {}

std::string Vec2InQuadInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Vec2InQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::vec2 Vec2InQuadInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inQuadInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inQuadInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InQuadInterpolationFunction::create() const {
    return new Vec2InQuadInterpolationFunction();
}

Vec2InCubicInterpolationFunction::Vec2InCubicInterpolationFunction() {}

std::string Vec2InCubicInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Vec2InCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::vec2 Vec2InCubicInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inCubicInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inCubicInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InCubicInterpolationFunction::create() const {
    return new Vec2InCubicInterpolationFunction();
}

Vec2InQuartInterpolationFunction::Vec2InQuartInterpolationFunction() {}

std::string Vec2InQuartInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Vec2InQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::vec2 Vec2InQuartInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inQuartInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inQuartInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InQuartInterpolationFunction::create() const {
    return new Vec2InQuartInterpolationFunction();
}

Vec2InQuintInterpolationFunction::Vec2InQuintInterpolationFunction() {}

std::string Vec2InQuintInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Vec2InQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::vec2 Vec2InQuintInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inQuintInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inQuintInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InQuintInterpolationFunction::create() const {
    return new Vec2InQuintInterpolationFunction();
}

Vec2InSineInterpolationFunction::Vec2InSineInterpolationFunction() {}

std::string Vec2InSineInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Vec2InSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::vec2 Vec2InSineInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inSineInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inSineInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InSineInterpolationFunction::create() const {
    return new Vec2InSineInterpolationFunction();
}

Vec2InExponentInterpolationFunction::Vec2InExponentInterpolationFunction() {}

std::string Vec2InExponentInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Vec2InExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::vec2 Vec2InExponentInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inExponentInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inExponentInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InExponentInterpolationFunction::create() const {
    return new Vec2InExponentInterpolationFunction();
}

Vec2InCircInterpolationFunction::Vec2InCircInterpolationFunction() {}

std::string Vec2InCircInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string Vec2InCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::vec2 Vec2InCircInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inCircInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inCircInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InCircInterpolationFunction::create() const {
    return new Vec2InCircInterpolationFunction();
}

Vec2OutQuadInterpolationFunction::Vec2OutQuadInterpolationFunction() {}

std::string Vec2OutQuadInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Vec2OutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::vec2 Vec2OutQuadInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outQuadInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outQuadInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutQuadInterpolationFunction::create() const {
    return new Vec2OutQuadInterpolationFunction();
}

Vec2OutCubicInterpolationFunction::Vec2OutCubicInterpolationFunction() {}

std::string Vec2OutCubicInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Vec2OutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::vec2 Vec2OutCubicInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outCubicInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outCubicInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutCubicInterpolationFunction::create() const {
    return new Vec2OutCubicInterpolationFunction();
}

Vec2OutQuartInterpolationFunction::Vec2OutQuartInterpolationFunction() {}

std::string Vec2OutQuartInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Vec2OutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::vec2 Vec2OutQuartInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outQuartInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outQuartInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutQuartInterpolationFunction::create() const {
    return new Vec2OutQuartInterpolationFunction();
}

Vec2OutQuintInterpolationFunction::Vec2OutQuintInterpolationFunction() {}

std::string Vec2OutQuintInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Vec2OutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::vec2 Vec2OutQuintInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outQuintInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outQuintInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutQuintInterpolationFunction::create() const {
    return new Vec2OutQuintInterpolationFunction();
}

Vec2OutSineInterpolationFunction::Vec2OutSineInterpolationFunction() {}

std::string Vec2OutSineInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Vec2OutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::vec2 Vec2OutSineInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outSineInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outSineInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutSineInterpolationFunction::create() const {
    return new Vec2OutSineInterpolationFunction();
}

Vec2OutExponentInterpolationFunction::Vec2OutExponentInterpolationFunction() {}

std::string Vec2OutExponentInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Vec2OutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::vec2 Vec2OutExponentInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outExponentInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outExponentInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutExponentInterpolationFunction::create() const {
    return new Vec2OutExponentInterpolationFunction();
}

Vec2OutCircInterpolationFunction::Vec2OutCircInterpolationFunction() {}

std::string Vec2OutCircInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string Vec2OutCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::vec2 Vec2OutCircInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outCircInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outCircInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutCircInterpolationFunction::create() const {
    return new Vec2OutCircInterpolationFunction();
}

Vec2InOutQuadInterpolationFunction::Vec2InOutQuadInterpolationFunction() {}

std::string Vec2InOutQuadInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Vec2InOutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::vec2 Vec2InOutQuadInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inOutQuadInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inOutQuadInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InOutQuadInterpolationFunction::create() const {
    return new Vec2InOutQuadInterpolationFunction();
}

Vec2InOutCubicInterpolationFunction::Vec2InOutCubicInterpolationFunction() {}

std::string Vec2InOutCubicInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Vec2InOutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::vec2 Vec2InOutCubicInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inOutCubicInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inOutCubicInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InOutCubicInterpolationFunction::create() const {
    return new Vec2InOutCubicInterpolationFunction();
}

Vec2InOutQuartInterpolationFunction::Vec2InOutQuartInterpolationFunction() {}

std::string Vec2InOutQuartInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Vec2InOutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::vec2 Vec2InOutQuartInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inOutQuartInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inOutQuartInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InOutQuartInterpolationFunction::create() const {
    return new Vec2InOutQuartInterpolationFunction();
}

Vec2InOutQuintInterpolationFunction::Vec2InOutQuintInterpolationFunction() {}

std::string Vec2InOutQuintInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Vec2InOutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::vec2 Vec2InOutQuintInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inOutQuintInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inOutQuintInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InOutQuintInterpolationFunction::create() const {
    return new Vec2InOutQuintInterpolationFunction();
}

Vec2InOutSineInterpolationFunction::Vec2InOutSineInterpolationFunction() {}

std::string Vec2InOutSineInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Vec2InOutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::vec2 Vec2InOutSineInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inOutSineInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inOutSineInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InOutSineInterpolationFunction::create() const {
    return new Vec2InOutSineInterpolationFunction();
}

Vec2InOutExponentInterpolationFunction::Vec2InOutExponentInterpolationFunction() {}

std::string Vec2InOutExponentInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Vec2InOutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::vec2 Vec2InOutExponentInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inOutExponentInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inOutExponentInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InOutExponentInterpolationFunction::create() const {
    return new Vec2InOutExponentInterpolationFunction();
}

Vec2InOutCircInterpolationFunction::Vec2InOutCircInterpolationFunction() {}

std::string Vec2InOutCircInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string Vec2InOutCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::vec2 Vec2InOutCircInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::inOutCircInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::inOutCircInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2InOutCircInterpolationFunction::create() const {
    return new Vec2InOutCircInterpolationFunction();
}

Vec2OutInQuadInterpolationFunction::Vec2OutInQuadInterpolationFunction() {}

std::string Vec2OutInQuadInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Vec2OutInQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::vec2 Vec2OutInQuadInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outInQuadInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outInQuadInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutInQuadInterpolationFunction::create() const {
    return new Vec2OutInQuadInterpolationFunction();
}

Vec2OutInCubicInterpolationFunction::Vec2OutInCubicInterpolationFunction() {}

std::string Vec2OutInCubicInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Vec2OutInCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::vec2 Vec2OutInCubicInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outInCubicInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outInCubicInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutInCubicInterpolationFunction::create() const {
    return new Vec2OutInCubicInterpolationFunction();
}

Vec2OutInQuartInterpolationFunction::Vec2OutInQuartInterpolationFunction() {}

std::string Vec2OutInQuartInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Vec2OutInQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::vec2 Vec2OutInQuartInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outInQuartInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outInQuartInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutInQuartInterpolationFunction::create() const {
    return new Vec2OutInQuartInterpolationFunction();
}

Vec2OutInQuintInterpolationFunction::Vec2OutInQuintInterpolationFunction() {}

std::string Vec2OutInQuintInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Vec2OutInQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::vec2 Vec2OutInQuintInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outInQuintInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outInQuintInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutInQuintInterpolationFunction::create() const {
    return new Vec2OutInQuintInterpolationFunction();
}

Vec2OutInSineInterpolationFunction::Vec2OutInSineInterpolationFunction() {}

std::string Vec2OutInSineInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Vec2OutInSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::vec2 Vec2OutInSineInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outInSineInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outInSineInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutInSineInterpolationFunction::create() const {
    return new Vec2OutInSineInterpolationFunction();
}

Vec2OutInExponentInterpolationFunction::Vec2OutInExponentInterpolationFunction() {}

std::string Vec2OutInExponentInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Vec2OutInExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::vec2 Vec2OutInExponentInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outInExponentInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outInExponentInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutInExponentInterpolationFunction::create() const {
    return new Vec2OutInExponentInterpolationFunction();
}

Vec2OutInCircInterpolationFunction::Vec2OutInCircInterpolationFunction() {}

std::string Vec2OutInCircInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string Vec2OutInCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::vec2 Vec2OutInCircInterpolationFunction::interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const {
    float compX = BasicFloatInterpolation::outInCircInterpolation(startvalue.x, endvalue.x, time);
    float compY = BasicFloatInterpolation::outInCircInterpolation(startvalue.y, endvalue.y, time);
    return tgt::vec2(compX, compY);
}

InterpolationFunction<tgt::vec2>* Vec2OutInCircInterpolationFunction::create() const {
    return new Vec2OutInCircInterpolationFunction();
}

} // namespace voreen
