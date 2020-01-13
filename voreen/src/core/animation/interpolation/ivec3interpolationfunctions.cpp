/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/core/animation/interpolation/basicintinterpolation.h"
#include "voreen/core/animation/interpolation/ivec3interpolationfunctions.h"

namespace voreen {

IVec3StartInterpolationFunction::IVec3StartInterpolationFunction() {}

std::string IVec3StartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string IVec3StartInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::ivec3 IVec3StartInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    if (time < 1.f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<tgt::ivec3>* IVec3StartInterpolationFunction::create() const {
    return new IVec3StartInterpolationFunction();
}

IVec3EndInterpolationFunction::IVec3EndInterpolationFunction() {}

std::string IVec3EndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}

std::string IVec3EndInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::ivec3 IVec3EndInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    if (time > 0.f)
        return endvalue;
    else
        return startvalue;
}

InterpolationFunction<tgt::ivec3>* IVec3EndInterpolationFunction::create() const {
    return new IVec3EndInterpolationFunction();
}

IVec3StartEndInterpolationFunction::IVec3StartEndInterpolationFunction() {}

std::string IVec3StartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string IVec3StartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

tgt::ivec3 IVec3StartEndInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    if (time < 0.5f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<tgt::ivec3>* IVec3StartEndInterpolationFunction::create() const {
    return new IVec3StartEndInterpolationFunction();
}

IVec3LinearInterpolationFunction::IVec3LinearInterpolationFunction() {}

std::string IVec3LinearInterpolationFunction::getGuiName() const {
    return "linear interpolation";
}

std::string IVec3LinearInterpolationFunction::getCategory() const {
    return "linear";
}

tgt::ivec3 IVec3LinearInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::linearInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::linearInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::linearInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3LinearInterpolationFunction::create() const {
    return new IVec3LinearInterpolationFunction();
}

IVec3InQuadInterpolationFunction::IVec3InQuadInterpolationFunction() {}

std::string IVec3InQuadInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IVec3InQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::ivec3 IVec3InQuadInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inQuadInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inQuadInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inQuadInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InQuadInterpolationFunction::create() const {
    return new IVec3InQuadInterpolationFunction();
}

IVec3InCubicInterpolationFunction::IVec3InCubicInterpolationFunction() {}

std::string IVec3InCubicInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IVec3InCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::ivec3 IVec3InCubicInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inCubicInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inCubicInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inCubicInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InCubicInterpolationFunction::create() const {
    return new IVec3InCubicInterpolationFunction();
}

IVec3InQuartInterpolationFunction::IVec3InQuartInterpolationFunction() {}

std::string IVec3InQuartInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IVec3InQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::ivec3 IVec3InQuartInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inQuartInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inQuartInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inQuartInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InQuartInterpolationFunction::create() const {
    return new IVec3InQuartInterpolationFunction();
}

IVec3InQuintInterpolationFunction::IVec3InQuintInterpolationFunction() {}

std::string IVec3InQuintInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IVec3InQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::ivec3 IVec3InQuintInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inQuintInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inQuintInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inQuintInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InQuintInterpolationFunction::create() const {
    return new IVec3InQuintInterpolationFunction();
}

IVec3InSineInterpolationFunction::IVec3InSineInterpolationFunction() {}

std::string IVec3InSineInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IVec3InSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::ivec3 IVec3InSineInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inSineInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inSineInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inSineInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InSineInterpolationFunction::create() const {
    return new IVec3InSineInterpolationFunction();
}

IVec3InExponentInterpolationFunction::IVec3InExponentInterpolationFunction() {}

std::string IVec3InExponentInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IVec3InExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::ivec3 IVec3InExponentInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inExponentInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inExponentInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inExponentInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InExponentInterpolationFunction::create() const {
    return new IVec3InExponentInterpolationFunction();
}

IVec3InCircInterpolationFunction::IVec3InCircInterpolationFunction() {}

std::string IVec3InCircInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IVec3InCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::ivec3 IVec3InCircInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inCircInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inCircInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inCircInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InCircInterpolationFunction::create() const {
    return new IVec3InCircInterpolationFunction();
}

IVec3OutQuadInterpolationFunction::IVec3OutQuadInterpolationFunction() {}

std::string IVec3OutQuadInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IVec3OutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::ivec3 IVec3OutQuadInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outQuadInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outQuadInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outQuadInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutQuadInterpolationFunction::create() const {
    return new IVec3OutQuadInterpolationFunction();
}

IVec3OutCubicInterpolationFunction::IVec3OutCubicInterpolationFunction() {}

std::string IVec3OutCubicInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IVec3OutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::ivec3 IVec3OutCubicInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outCubicInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outCubicInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outCubicInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutCubicInterpolationFunction::create() const {
    return new IVec3OutCubicInterpolationFunction();
}

IVec3OutQuartInterpolationFunction::IVec3OutQuartInterpolationFunction() {}

std::string IVec3OutQuartInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IVec3OutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::ivec3 IVec3OutQuartInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outQuartInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outQuartInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outQuartInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutQuartInterpolationFunction::create() const {
    return new IVec3OutQuartInterpolationFunction();
}

IVec3OutQuintInterpolationFunction::IVec3OutQuintInterpolationFunction() {}

std::string IVec3OutQuintInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IVec3OutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::ivec3 IVec3OutQuintInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outQuintInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outQuintInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outQuintInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutQuintInterpolationFunction::create() const {
    return new IVec3OutQuintInterpolationFunction();
}

IVec3OutSineInterpolationFunction::IVec3OutSineInterpolationFunction() {}

std::string IVec3OutSineInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IVec3OutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::ivec3 IVec3OutSineInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outSineInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outSineInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outSineInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutSineInterpolationFunction::create() const {
    return new IVec3OutSineInterpolationFunction();
}

IVec3OutExponentInterpolationFunction::IVec3OutExponentInterpolationFunction() {}

std::string IVec3OutExponentInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IVec3OutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::ivec3 IVec3OutExponentInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outExponentInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outExponentInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outExponentInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutExponentInterpolationFunction::create() const {
    return new IVec3OutExponentInterpolationFunction();
}

IVec3OutCircInterpolationFunction::IVec3OutCircInterpolationFunction() {}

std::string IVec3OutCircInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IVec3OutCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::ivec3 IVec3OutCircInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outCircInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outCircInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outCircInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutCircInterpolationFunction::create() const {
    return new IVec3OutCircInterpolationFunction();
}

IVec3InOutQuadInterpolationFunction::IVec3InOutQuadInterpolationFunction() {}

std::string IVec3InOutQuadInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IVec3InOutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::ivec3 IVec3InOutQuadInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inOutQuadInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inOutQuadInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inOutQuadInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InOutQuadInterpolationFunction::create() const {
    return new IVec3InOutQuadInterpolationFunction();
}

IVec3InOutCubicInterpolationFunction::IVec3InOutCubicInterpolationFunction() {}

std::string IVec3InOutCubicInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IVec3InOutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::ivec3 IVec3InOutCubicInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inOutCubicInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inOutCubicInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inOutCubicInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InOutCubicInterpolationFunction::create() const {
    return new IVec3InOutCubicInterpolationFunction();
}

IVec3InOutQuartInterpolationFunction::IVec3InOutQuartInterpolationFunction() {}

std::string IVec3InOutQuartInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IVec3InOutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::ivec3 IVec3InOutQuartInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inOutQuartInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inOutQuartInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inOutQuartInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InOutQuartInterpolationFunction::create() const {
    return new IVec3InOutQuartInterpolationFunction();
}

IVec3InOutQuintInterpolationFunction::IVec3InOutQuintInterpolationFunction() {}

std::string IVec3InOutQuintInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IVec3InOutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::ivec3 IVec3InOutQuintInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inOutQuintInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inOutQuintInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inOutQuintInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InOutQuintInterpolationFunction::create() const {
    return new IVec3InOutQuintInterpolationFunction();
}

IVec3InOutSineInterpolationFunction::IVec3InOutSineInterpolationFunction() {}

std::string IVec3InOutSineInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IVec3InOutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::ivec3 IVec3InOutSineInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inOutSineInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inOutSineInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inOutSineInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InOutSineInterpolationFunction::create() const {
    return new IVec3InOutSineInterpolationFunction();
}

IVec3InOutExponentInterpolationFunction::IVec3InOutExponentInterpolationFunction() {}

std::string IVec3InOutExponentInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IVec3InOutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::ivec3 IVec3InOutExponentInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inOutExponentInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inOutExponentInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inOutExponentInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InOutExponentInterpolationFunction::create() const {
    return new IVec3InOutExponentInterpolationFunction();
}

IVec3InOutCircInterpolationFunction::IVec3InOutCircInterpolationFunction() {}

std::string IVec3InOutCircInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IVec3InOutCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::ivec3 IVec3InOutCircInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::inOutCircInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::inOutCircInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::inOutCircInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3InOutCircInterpolationFunction::create() const {
    return new IVec3InOutCircInterpolationFunction();
}

IVec3OutInQuadInterpolationFunction::IVec3OutInQuadInterpolationFunction() {}

std::string IVec3OutInQuadInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IVec3OutInQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

tgt::ivec3 IVec3OutInQuadInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outInQuadInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outInQuadInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outInQuadInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutInQuadInterpolationFunction::create() const {
    return new IVec3OutInQuadInterpolationFunction();
}

IVec3OutInCubicInterpolationFunction::IVec3OutInCubicInterpolationFunction() {}

std::string IVec3OutInCubicInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IVec3OutInCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

tgt::ivec3 IVec3OutInCubicInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outInCubicInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outInCubicInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outInCubicInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutInCubicInterpolationFunction::create() const {
    return new IVec3OutInCubicInterpolationFunction();
}

IVec3OutInQuartInterpolationFunction::IVec3OutInQuartInterpolationFunction() {}

std::string IVec3OutInQuartInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IVec3OutInQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

tgt::ivec3 IVec3OutInQuartInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outInQuartInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outInQuartInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outInQuartInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutInQuartInterpolationFunction::create() const {
    return new IVec3OutInQuartInterpolationFunction();
}

IVec3OutInQuintInterpolationFunction::IVec3OutInQuintInterpolationFunction() {}

std::string IVec3OutInQuintInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IVec3OutInQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

tgt::ivec3 IVec3OutInQuintInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outInQuintInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outInQuintInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outInQuintInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutInQuintInterpolationFunction::create() const {
    return new IVec3OutInQuintInterpolationFunction();
}

IVec3OutInSineInterpolationFunction::IVec3OutInSineInterpolationFunction() {}

std::string IVec3OutInSineInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IVec3OutInSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

tgt::ivec3 IVec3OutInSineInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outInSineInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outInSineInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outInSineInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutInSineInterpolationFunction::create() const {
    return new IVec3OutInSineInterpolationFunction();
}

IVec3OutInExponentInterpolationFunction::IVec3OutInExponentInterpolationFunction() {}

std::string IVec3OutInExponentInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IVec3OutInExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

tgt::ivec3 IVec3OutInExponentInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outInExponentInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outInExponentInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outInExponentInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutInExponentInterpolationFunction::create() const {
    return new IVec3OutInExponentInterpolationFunction();
}

IVec3OutInCircInterpolationFunction::IVec3OutInCircInterpolationFunction() {}

std::string IVec3OutInCircInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IVec3OutInCircInterpolationFunction::getCategory() const {
    return "circular";
}

tgt::ivec3 IVec3OutInCircInterpolationFunction::interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const {
    int compX = BasicIntInterpolation::outInCircInterpolation(static_cast<float>(startvalue.x), static_cast<float>(endvalue.x), time);
    int compY = BasicIntInterpolation::outInCircInterpolation(static_cast<float>(startvalue.y), static_cast<float>(endvalue.y), time);
    int compZ = BasicIntInterpolation::outInCircInterpolation(static_cast<float>(startvalue.z), static_cast<float>(endvalue.z), time);
    return tgt::ivec3(compX, compY, compZ);
}

InterpolationFunction<tgt::ivec3>* IVec3OutInCircInterpolationFunction::create() const {
    return new IVec3OutInCircInterpolationFunction();
}

} // namespace voreen
