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

#include "voreen/core/animation/interpolation/basicintinterpolation.h"
#include "voreen/core/animation/interpolation/intinterpolationfunctions.h"
#include "voreen/core/animation/propertykeyvalue.h"

namespace voreen {

IntStartInterpolationFunction::IntStartInterpolationFunction() {}

std::string IntStartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string IntStartInterpolationFunction::getCategory() const {
    return "boolean";
}

int IntStartInterpolationFunction::interpolate(int startvalue, int endvalue, float time) const {
    if (time < 1.f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<int>* IntStartInterpolationFunction::create() const {
    return new IntStartInterpolationFunction();
}

IntEndInterpolationFunction::IntEndInterpolationFunction() {}

std::string IntEndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}

std::string IntEndInterpolationFunction::getCategory() const {
    return "boolean";
}

int IntEndInterpolationFunction::interpolate(int startvalue, int endvalue, float time) const {
    if (time > 0.f)
        return endvalue;
    else
        return startvalue;
}
InterpolationFunction<int>* IntEndInterpolationFunction::create() const {
    return new IntEndInterpolationFunction();
}

IntStartEndInterpolationFunction::IntStartEndInterpolationFunction() {}

std::string IntStartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string IntStartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

int IntStartEndInterpolationFunction::interpolate(int startvalue, int endvalue, float time) const {
    if (time < 0.5f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<int>* IntStartEndInterpolationFunction::create() const {
    return new IntStartEndInterpolationFunction();
}

IntLinearInterpolationFunction::IntLinearInterpolationFunction() {}

std::string IntLinearInterpolationFunction::getGuiName() const {
    return "linear interpolation";
}

std::string IntLinearInterpolationFunction::getCategory() const {
    return "linear";
}

int IntLinearInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::linearInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntLinearInterpolationFunction::create() const {
    return new IntLinearInterpolationFunction();
}

IntCatmullRomInterpolationFunction::IntCatmullRomInterpolationFunction() {}

std::string IntCatmullRomInterpolationFunction::getGuiName() const {
        return "Catmull-Rom spline interpolation";
}

std::string IntCatmullRomInterpolationFunction::getCategory() const {
        return "spline";
}

int IntCatmullRomInterpolationFunction::interpolate(std::vector<PropertyKeyValue<int>*> controlpoints, float time) const {
    std::vector<PropertyKeyValue<int>*>::iterator it;
    it = controlpoints.begin();
    while ((*it)->getTime() < time)
        it++;

    float time2 = (*it)->getTime();
    int p2 = (*it)->getValue();
    it--;
    float time1 = (*it)->getTime();
    int p1 = (*it)->getValue();
    int m1;
    if (it!= controlpoints.begin()) {
        it--;
        int p0 = static_cast<int>((*it)->getValue());
        m1 = static_cast<int>(0.5f*(p2 - p0));
        it++;
    } else
        m1 = p2 - p1;

    it++;
    it++;
    int m2;
    if (it!= controlpoints.end()) {
        int p3 = (*it)->getValue();
        m2 = static_cast<int>(0.5f*(p3 - p1));
    } else
        m2 = p2 - p1;

    float t = (time - time1)/(time2 - time1);
    const float h00 = (1 + 2*t)*(1 - t)*(1 - t);
    const float h10 = t*(1 - t)*(1 - t);
    const float h01 = t*t*(3-2*t);
    const float h11 = t*t*(t - 1);
    return static_cast<int>(h00*p1 + h10*m1 + h01*p2 + h11*m2);
}

MultiPointInterpolationFunction<int>* IntCatmullRomInterpolationFunction::create() const {
        return new IntCatmullRomInterpolationFunction();
}

IntInQuadInterpolationFunction::IntInQuadInterpolationFunction() {}

std::string IntInQuadInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IntInQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

int IntInQuadInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inQuadInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInQuadInterpolationFunction::create() const {
    return new IntInQuadInterpolationFunction();
}

IntInCubicInterpolationFunction::IntInCubicInterpolationFunction() {}

std::string IntInCubicInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IntInCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

int IntInCubicInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inCubicInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInCubicInterpolationFunction::create() const {
    return new IntInCubicInterpolationFunction();
}

IntInQuartInterpolationFunction::IntInQuartInterpolationFunction() {}

std::string IntInQuartInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IntInQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

int IntInQuartInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inQuartInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInQuartInterpolationFunction::create() const {
    return new IntInQuartInterpolationFunction();
}

IntInQuintInterpolationFunction::IntInQuintInterpolationFunction() {}

std::string IntInQuintInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IntInQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

int IntInQuintInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inQuintInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInQuintInterpolationFunction::create() const {
    return new IntInQuintInterpolationFunction();
}

IntInSineInterpolationFunction::IntInSineInterpolationFunction() {}

std::string IntInSineInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IntInSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

int IntInSineInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inSineInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInSineInterpolationFunction::create() const {
    return new IntInSineInterpolationFunction();
}

IntInExponentInterpolationFunction::IntInExponentInterpolationFunction() {}

std::string IntInExponentInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IntInExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

int IntInExponentInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inExponentInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInExponentInterpolationFunction::create() const {
    return new IntInExponentInterpolationFunction();
}

IntInCircInterpolationFunction::IntInCircInterpolationFunction() {}

std::string IntInCircInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string IntInCircInterpolationFunction::getCategory() const {
    return "circular";
}

int IntInCircInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inCircInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInCircInterpolationFunction::create() const {
    return new IntInCircInterpolationFunction();
}

IntOutQuadInterpolationFunction::IntOutQuadInterpolationFunction() {}

std::string IntOutQuadInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IntOutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

int IntOutQuadInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outQuadInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutQuadInterpolationFunction::create() const {
    return new IntOutQuadInterpolationFunction();
}

IntOutCubicInterpolationFunction::IntOutCubicInterpolationFunction() {}

std::string IntOutCubicInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IntOutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

int IntOutCubicInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outCubicInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutCubicInterpolationFunction::create() const {
    return new IntOutCubicInterpolationFunction();
}

IntOutQuartInterpolationFunction::IntOutQuartInterpolationFunction() {}

std::string IntOutQuartInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IntOutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

int IntOutQuartInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outQuartInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutQuartInterpolationFunction::create() const {
    return new IntOutQuartInterpolationFunction();
}

IntOutQuintInterpolationFunction::IntOutQuintInterpolationFunction() {}

std::string IntOutQuintInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IntOutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

int IntOutQuintInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outQuintInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutQuintInterpolationFunction::create() const {
    return new IntOutQuintInterpolationFunction();
}

IntOutSineInterpolationFunction::IntOutSineInterpolationFunction() {}

std::string IntOutSineInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IntOutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

int IntOutSineInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outSineInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutSineInterpolationFunction::create() const {
    return new IntOutSineInterpolationFunction();
}

IntOutExponentInterpolationFunction::IntOutExponentInterpolationFunction() {}

std::string IntOutExponentInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IntOutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

int IntOutExponentInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outExponentInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutExponentInterpolationFunction::create() const {
    return new IntOutExponentInterpolationFunction();
}

IntOutCircInterpolationFunction::IntOutCircInterpolationFunction() {}

std::string IntOutCircInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string IntOutCircInterpolationFunction::getCategory() const {
    return "circular";
}

int IntOutCircInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outCircInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutCircInterpolationFunction::create() const {
    return new IntOutCircInterpolationFunction();
}

IntInOutQuadInterpolationFunction::IntInOutQuadInterpolationFunction() {}

std::string IntInOutQuadInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IntInOutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

int IntInOutQuadInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inOutQuadInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInOutQuadInterpolationFunction::create() const {
    return new IntInOutQuadInterpolationFunction();
}

IntInOutCubicInterpolationFunction::IntInOutCubicInterpolationFunction() {}

std::string IntInOutCubicInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IntInOutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

int IntInOutCubicInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inOutCubicInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInOutCubicInterpolationFunction::create() const {
    return new IntInOutCubicInterpolationFunction();
}

IntInOutQuartInterpolationFunction::IntInOutQuartInterpolationFunction() {}

std::string IntInOutQuartInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IntInOutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

int IntInOutQuartInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inOutQuartInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInOutQuartInterpolationFunction::create() const {
    return new IntInOutQuartInterpolationFunction();
}

IntInOutQuintInterpolationFunction::IntInOutQuintInterpolationFunction() {}

std::string IntInOutQuintInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IntInOutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

int IntInOutQuintInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inOutQuintInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInOutQuintInterpolationFunction::create() const {
    return new IntInOutQuintInterpolationFunction();
}

IntInOutSineInterpolationFunction::IntInOutSineInterpolationFunction() {}

std::string IntInOutSineInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IntInOutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

int IntInOutSineInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inOutSineInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInOutSineInterpolationFunction::create() const {
    return new IntInOutSineInterpolationFunction();
}

IntInOutExponentInterpolationFunction::IntInOutExponentInterpolationFunction() {}

std::string IntInOutExponentInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IntInOutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

int IntInOutExponentInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inOutExponentInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInOutExponentInterpolationFunction::create() const {
    return new IntInOutExponentInterpolationFunction();
}

IntInOutCircInterpolationFunction::IntInOutCircInterpolationFunction() {}

std::string IntInOutCircInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string IntInOutCircInterpolationFunction::getCategory() const {
    return "circular";
}

int IntInOutCircInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::inOutCircInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntInOutCircInterpolationFunction::create() const {
    return new IntInOutCircInterpolationFunction();
}

IntOutInQuadInterpolationFunction::IntOutInQuadInterpolationFunction() {}

std::string IntOutInQuadInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IntOutInQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

int IntOutInQuadInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outInQuadInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutInQuadInterpolationFunction::create() const {
    return new IntOutInQuadInterpolationFunction();
}

IntOutInCubicInterpolationFunction::IntOutInCubicInterpolationFunction() {}

std::string IntOutInCubicInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IntOutInCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

int IntOutInCubicInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outInCubicInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutInCubicInterpolationFunction::create() const {
    return new IntOutInCubicInterpolationFunction();
}

IntOutInQuartInterpolationFunction::IntOutInQuartInterpolationFunction() {}

std::string IntOutInQuartInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IntOutInQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

int IntOutInQuartInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outInQuartInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutInQuartInterpolationFunction::create() const {
    return new IntOutInQuartInterpolationFunction();
}

IntOutInQuintInterpolationFunction::IntOutInQuintInterpolationFunction() {}

std::string IntOutInQuintInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IntOutInQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

int IntOutInQuintInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outInQuintInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutInQuintInterpolationFunction::create() const {
    return new IntOutInQuintInterpolationFunction();
}

IntOutInSineInterpolationFunction::IntOutInSineInterpolationFunction() {}

std::string IntOutInSineInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IntOutInSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

int IntOutInSineInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outInSineInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutInSineInterpolationFunction::create() const {
    return new IntOutInSineInterpolationFunction();
}

IntOutInExponentInterpolationFunction::IntOutInExponentInterpolationFunction() {}

std::string IntOutInExponentInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IntOutInExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

int IntOutInExponentInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outInExponentInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutInExponentInterpolationFunction::create() const {
    return new IntOutInExponentInterpolationFunction();
}

IntOutInCircInterpolationFunction::IntOutInCircInterpolationFunction() {}

std::string IntOutInCircInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string IntOutInCircInterpolationFunction::getCategory() const {
    return "circular";
}

int IntOutInCircInterpolationFunction::interpolate(int startvalue,int endvalue, float time) const {
    return BasicIntInterpolation::outInCircInterpolation(static_cast<float>(startvalue), static_cast<float>(endvalue), time);
}

InterpolationFunction<int>* IntOutInCircInterpolationFunction::create() const {
    return new IntOutInCircInterpolationFunction();
}

} // namespace voreen
