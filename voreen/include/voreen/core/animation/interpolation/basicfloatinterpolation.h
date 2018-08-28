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

#ifndef VRN_BASICFLOATINTERPOLATION_H
#define VRN_BASICFLOATINTERPOLATION_H

#include "voreen/core/voreencoreapi.h"

namespace voreen {

/**
 * This class offers a set of basic interpolation functions for float-values.
 */
class VRN_CORE_API BasicFloatInterpolation {
public:
    BasicFloatInterpolation();

    /**
    * interpolates linear between the startvalue and the endvalue
    * timeinterval 0-1
    */
    static float linearInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function x^2
     * timeinterval 0-1
     */

    static float inQuadInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function x^3
     * timeinterval 0-1
     */
    static float inCubicInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function x^4
     * timeinterval 0-1
     */
    static float inQuartInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function x^5
     * timeinterval 0-1
     */
    static float inQuintInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function sin(x/2*PI)
     * timeinterval 0-1
     */
    static float inSineInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with an exponential function 2^(10*(x-1))-0.0001
     * timeinterval 0-1
     */
    static float inExponentInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function 1-(sqrt(1-x^2))
     * timeinterval 0-1
     */
    static float inCircInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuadInterpolation, but easing out
     * timeinterval 0-1
     */
    static float outQuadInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCubicInterpolation, but easing out
     * timeinterval 0-1
     */
    static float outCubicInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuartInterpolation, but easing out
     * timeinterval 0-1
     */
    static float outQuartInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuintInterpolation, but easing out
     * timeinterval 0-1
     */
    static float outQuintInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inSinusInterpolation, but easing out
     * timeinterval 0-1
     */
    static float outSineInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inExponentInterpolation, but easing out
     * timeinterval 0-1
     */
    static float outExponentInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCircInterpolation, but easing out
     * timeinterval 0-1
     */
    static float outCircInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuadInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static float inOutQuadInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCubicInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static float inOutCubicInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuartInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static float inOutQuartInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuintInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static float inOutQuintInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inSinusInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static float inOutSineInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inExponentInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static float inOutExponentInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCircInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static float inOutCircInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuadInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static float outInQuadInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCubicInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static float outInCubicInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuartInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static float outInQuartInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuintInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static float outInQuintInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inSinusInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static float outInSineInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inExponentInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static float outInExponentInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCircInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static float outInCircInterpolation(float startvalue, float endvalue, float time);
};

} // namespace voreen
#endif
