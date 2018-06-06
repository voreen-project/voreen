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
#include "voreen/core/animation/interpolation/floatinterpolationfunctions.h"
#include "voreen/core/animation/propertykeyvalue.h"

namespace voreen {

FloatStartInterpolationFunction::FloatStartInterpolationFunction() {}

std::string FloatStartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string FloatStartInterpolationFunction::getCategory() const {
    return "boolean";
}

float FloatStartInterpolationFunction::interpolate(float startvalue, float endvalue, float time) const {
    if (time<1)
        return startvalue;
    else
        return endvalue;
}
InterpolationFunction<float>* FloatStartInterpolationFunction::create() const {
    return new FloatStartInterpolationFunction();
}

FloatEndInterpolationFunction::FloatEndInterpolationFunction() {}

std::string FloatEndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}

std::string FloatEndInterpolationFunction::getCategory() const {
    return "boolean";
}

float FloatEndInterpolationFunction::interpolate(float startvalue, float endvalue, float time) const {
    if (time > 0.f)
        return endvalue;
    else
        return startvalue;
}

InterpolationFunction<float>* FloatEndInterpolationFunction::create() const {
    return new FloatEndInterpolationFunction();
}

FloatStartEndInterpolationFunction::FloatStartEndInterpolationFunction() {}

std::string FloatStartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string FloatStartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

float FloatStartEndInterpolationFunction::interpolate(float startvalue, float endvalue, float time) const {
    if (time<0.5){
        return startvalue;
    }
    else{
        return endvalue;
    }
}

InterpolationFunction<float>* FloatStartEndInterpolationFunction::create() const {
    return new FloatStartEndInterpolationFunction();
}

FloatLinearInterpolationFunction::FloatLinearInterpolationFunction() {}

std::string FloatLinearInterpolationFunction::getGuiName() const {
    return "linear interpolation";
}

std::string FloatLinearInterpolationFunction::getCategory() const {
    return "linear";
}

float FloatLinearInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::linearInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatLinearInterpolationFunction::create() const {
    return new FloatLinearInterpolationFunction();
}

FloatCatmullRomInterpolationFunction::FloatCatmullRomInterpolationFunction() {}

std::string FloatCatmullRomInterpolationFunction::getGuiName() const {
        return "Catmull-Rom spline interpolation";
}

std::string FloatCatmullRomInterpolationFunction::getCategory() const {
        return "spline";
}

float FloatCatmullRomInterpolationFunction::interpolate(std::vector<PropertyKeyValue<float>*> controlpoints, float time) const {
    std::vector<PropertyKeyValue<float>*>::iterator it;
    it = controlpoints.begin();
    while ((*it)->getTime() < time)
        it++;

    float time2 = (*it)->getTime();
    float p2 = (*it)->getValue();
    it--;
    float time1 = (*it)->getTime();
    float p1 = (*it)->getValue();
    float m1;
    if (it!= controlpoints.begin()){
        it--;
        float p0 = (*it)->getValue();
        m1 = 0.5f*(p2 - p0);
        it++;
    } else
        m1 = p2 - p1;

    it++;
    it++;
    float m2;
    if (it!= controlpoints.end()) {
        float p3 = (*it)->getValue();
        m2 = 0.5f*(p3 - p1);
    } else
        m2 = p2 - p1;

    float t = (time - time1)/(time2 - time1);
    const float h00 = (1 + 2*t)*(1 - t)*(1 - t);
    const float h10 = t*(1 - t)*(1 - t);
    const float h01 = t*t*(3-2*t);
    const float h11 = t*t*(t - 1);
    return h00*p1 + h10*m1 + h01*p2 + h11*m2;
}

MultiPointInterpolationFunction<float>* FloatCatmullRomInterpolationFunction::create() const {
    return new FloatCatmullRomInterpolationFunction();
}

FloatInQuadInterpolationFunction::FloatInQuadInterpolationFunction() {}

std::string FloatInQuadInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string FloatInQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

float FloatInQuadInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inQuadInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInQuadInterpolationFunction::create() const {
    return new FloatInQuadInterpolationFunction();
}

FloatInCubicInterpolationFunction::FloatInCubicInterpolationFunction() {}

std::string FloatInCubicInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string FloatInCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

float FloatInCubicInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inCubicInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInCubicInterpolationFunction::create() const {
    return new FloatInCubicInterpolationFunction();
}

FloatInQuartInterpolationFunction::FloatInQuartInterpolationFunction() {}

std::string FloatInQuartInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string FloatInQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

float FloatInQuartInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inQuartInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInQuartInterpolationFunction::create() const {
    return new FloatInQuartInterpolationFunction();
}

FloatInQuintInterpolationFunction::FloatInQuintInterpolationFunction() {}

std::string FloatInQuintInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string FloatInQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

float FloatInQuintInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inQuintInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInQuintInterpolationFunction::create() const {
    return new FloatInQuintInterpolationFunction();
}

FloatInSineInterpolationFunction::FloatInSineInterpolationFunction() {}

std::string FloatInSineInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string FloatInSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

float FloatInSineInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inSineInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInSineInterpolationFunction::create() const {
    return new FloatInSineInterpolationFunction();
}

FloatInExponentInterpolationFunction::FloatInExponentInterpolationFunction() {}

std::string FloatInExponentInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string FloatInExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

float FloatInExponentInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inExponentInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInExponentInterpolationFunction::create() const {
    return new FloatInExponentInterpolationFunction();
}

FloatInCircInterpolationFunction::FloatInCircInterpolationFunction() {}

std::string FloatInCircInterpolationFunction::getGuiName() const {
    return "easing in";
}

std::string FloatInCircInterpolationFunction::getCategory() const {
    return "circular";
}

float FloatInCircInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inCircInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInCircInterpolationFunction::create() const {
    return new FloatInCircInterpolationFunction();
}

FloatOutQuadInterpolationFunction::FloatOutQuadInterpolationFunction() {}

std::string FloatOutQuadInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string FloatOutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

float FloatOutQuadInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outQuadInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutQuadInterpolationFunction::create() const {
    return new FloatOutQuadInterpolationFunction();
}

FloatOutCubicInterpolationFunction::FloatOutCubicInterpolationFunction() {}

std::string FloatOutCubicInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string FloatOutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

float FloatOutCubicInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outCubicInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutCubicInterpolationFunction::create() const {
    return new FloatOutCubicInterpolationFunction();
}

FloatOutQuartInterpolationFunction::FloatOutQuartInterpolationFunction() {}

std::string FloatOutQuartInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string FloatOutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

float FloatOutQuartInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outQuartInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutQuartInterpolationFunction::create() const {
    return new FloatOutQuartInterpolationFunction();
}

FloatOutQuintInterpolationFunction::FloatOutQuintInterpolationFunction() {}

std::string FloatOutQuintInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string FloatOutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

float FloatOutQuintInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outQuintInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutQuintInterpolationFunction::create() const {
    return new FloatOutQuintInterpolationFunction();
}

FloatOutSineInterpolationFunction::FloatOutSineInterpolationFunction() {}

std::string FloatOutSineInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string FloatOutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

float FloatOutSineInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outSineInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutSineInterpolationFunction::create() const {
    return new FloatOutSineInterpolationFunction();
}

FloatOutExponentInterpolationFunction::FloatOutExponentInterpolationFunction() {}

std::string FloatOutExponentInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string FloatOutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

float FloatOutExponentInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outExponentInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutExponentInterpolationFunction::create() const {
    return new FloatOutExponentInterpolationFunction();
}

FloatOutCircInterpolationFunction::FloatOutCircInterpolationFunction() {}

std::string FloatOutCircInterpolationFunction::getGuiName() const {
    return "easing out";
}

std::string FloatOutCircInterpolationFunction::getCategory() const {
    return "circular";
}

float FloatOutCircInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outCircInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutCircInterpolationFunction::create() const {
    return new FloatOutCircInterpolationFunction();
}

FloatInOutQuadInterpolationFunction::FloatInOutQuadInterpolationFunction() {}

std::string FloatInOutQuadInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string FloatInOutQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

float FloatInOutQuadInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inOutQuadInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInOutQuadInterpolationFunction::create() const {
    return new FloatInOutQuadInterpolationFunction();
}

FloatInOutCubicInterpolationFunction::FloatInOutCubicInterpolationFunction() {}

std::string FloatInOutCubicInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string FloatInOutCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

float FloatInOutCubicInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inOutCubicInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInOutCubicInterpolationFunction::create() const {
    return new FloatInOutCubicInterpolationFunction();
}

FloatInOutQuartInterpolationFunction::FloatInOutQuartInterpolationFunction() {}

std::string FloatInOutQuartInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string FloatInOutQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

float FloatInOutQuartInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inOutQuartInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInOutQuartInterpolationFunction::create() const {
    return new FloatInOutQuartInterpolationFunction();
}

FloatInOutQuintInterpolationFunction::FloatInOutQuintInterpolationFunction() {}

std::string FloatInOutQuintInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string FloatInOutQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

float FloatInOutQuintInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inOutQuintInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInOutQuintInterpolationFunction::create() const {
    return new FloatInOutQuintInterpolationFunction();
}

FloatInOutSineInterpolationFunction::FloatInOutSineInterpolationFunction() {}

std::string FloatInOutSineInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string FloatInOutSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

float FloatInOutSineInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inOutSineInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInOutSineInterpolationFunction::create() const {
    return new FloatInOutSineInterpolationFunction();
}

FloatInOutExponentInterpolationFunction::FloatInOutExponentInterpolationFunction() {}

std::string FloatInOutExponentInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string FloatInOutExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

float FloatInOutExponentInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inOutExponentInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInOutExponentInterpolationFunction::create() const {
    return new FloatInOutExponentInterpolationFunction();
}

FloatInOutCircInterpolationFunction::FloatInOutCircInterpolationFunction() {}

std::string FloatInOutCircInterpolationFunction::getGuiName() const {
    return "first easing in, then easing out";
}

std::string FloatInOutCircInterpolationFunction::getCategory() const {
    return "circular";
}

float FloatInOutCircInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::inOutCircInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatInOutCircInterpolationFunction::create() const {
    return new FloatInOutCircInterpolationFunction();
}

FloatOutInQuadInterpolationFunction::FloatOutInQuadInterpolationFunction() {}

std::string FloatOutInQuadInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string FloatOutInQuadInterpolationFunction::getCategory() const {
    return "quadratic";
}

float FloatOutInQuadInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outInQuadInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutInQuadInterpolationFunction::create() const {
    return new FloatOutInQuadInterpolationFunction();
}

FloatOutInCubicInterpolationFunction::FloatOutInCubicInterpolationFunction() {}

std::string FloatOutInCubicInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string FloatOutInCubicInterpolationFunction::getCategory() const {
    return "cubic";
}

float FloatOutInCubicInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outInCubicInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutInCubicInterpolationFunction::create() const {
    return new FloatOutInCubicInterpolationFunction();
}

FloatOutInQuartInterpolationFunction::FloatOutInQuartInterpolationFunction() {}

std::string FloatOutInQuartInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string FloatOutInQuartInterpolationFunction::getCategory() const {
    return "quartetic";
}

float FloatOutInQuartInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outInQuartInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutInQuartInterpolationFunction::create() const {
    return new FloatOutInQuartInterpolationFunction();
}

FloatOutInQuintInterpolationFunction::FloatOutInQuintInterpolationFunction() {}

std::string FloatOutInQuintInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string FloatOutInQuintInterpolationFunction::getCategory() const {
    return "quintic";
}

float FloatOutInQuintInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outInQuintInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutInQuintInterpolationFunction::create() const {
    return new FloatOutInQuintInterpolationFunction();
}

FloatOutInSineInterpolationFunction::FloatOutInSineInterpolationFunction() {}

std::string FloatOutInSineInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string FloatOutInSineInterpolationFunction::getCategory() const {
    return "sineousidal";
}

float FloatOutInSineInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outInSineInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutInSineInterpolationFunction::create() const {
    return new FloatOutInSineInterpolationFunction();
}

FloatOutInExponentInterpolationFunction::FloatOutInExponentInterpolationFunction() {}

std::string FloatOutInExponentInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string FloatOutInExponentInterpolationFunction::getCategory() const {
    return "exponential";
}

float FloatOutInExponentInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outInExponentInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutInExponentInterpolationFunction::create() const {
    return new FloatOutInExponentInterpolationFunction();
}

FloatOutInCircInterpolationFunction::FloatOutInCircInterpolationFunction() {}

std::string FloatOutInCircInterpolationFunction::getGuiName() const {
    return "first easing out, then easing in";
}

std::string FloatOutInCircInterpolationFunction::getCategory() const {
    return "circular";
}

float FloatOutInCircInterpolationFunction::interpolate(float startvalue,float endvalue, float time) const {
    return BasicFloatInterpolation::outInCircInterpolation(startvalue,endvalue,time);
}

InterpolationFunction<float>* FloatOutInCircInterpolationFunction::create() const {
    return new FloatOutInCircInterpolationFunction();
}

} // namespace voreen
