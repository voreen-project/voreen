/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_BASICINTINTERPOLATION_H
#define VRN_BASICINTINTERPOLATION_H

namespace voreen {

/**
 * This class offers a set of basic interpolation functions for int-values.
 */
class BasicIntInterpolation {
public:
    BasicIntInterpolation();

    /**
     * interpolates linear between the startvalue and the endvalue
     * timeinterval 0-1
     */
    static int linearInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function x^2
     * timeinterval 0-1
     */
    static int inQuadInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function x^3
     * timeinterval 0-1
     */
    static int inCubicInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function x^4
     * timeinterval 0-1
     */
    static int inQuartInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function x^5
     * timeinterval 0-1
     */
    static int inQuintInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function sin(x/2*PI)
     * timeinterval 0-1
     */
    static int inSineInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with an exponentialfunction 2^(10*(x-1))-0.0001
     * timeinterval 0-1
     */
    static int inExponentInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue with function 1-(sqrt(1-x^2))
     * timeinterval 0-1
     */
    static int inCircInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuadInterpolation, but easing out
     * timeinterval 0-1
     */
    static int outQuadInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCubicInterpolation, but easing out
     * timeinterval 0-1
     */
    static int outCubicInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuartInterpolation, but easing out
     * timeinterval 0-1
     */
    static int outQuartInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuintInterpolation, but easing out
     * timeinterval 0-1
     */
    static int outQuintInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inSinusInterpolation, but easing out
     * timeinterval 0-1
     */
    static int outSineInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inExponentInterpolation, but easing out
     * timeinterval 0-1
     */
    static int outExponentInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCircInterpolation, but easing out
     * timeinterval 0-1
     */
    static int outCircInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuadInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static int inOutQuadInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCubicInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static int inOutCubicInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuartInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static int inOutQuartInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuintInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static int inOutQuintInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inSinusInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static int inOutSineInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inExponentInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static int inOutExponentInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCircInterpolation, but first easing in until x=0.5 then easing out
     * timeinterval 0-1
     */
    static int inOutCircInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuadInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static int outInQuadInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCubicInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static int outInCubicInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuartInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static int outInQuartInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inQuintInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static int outInQuintInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inSinusInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static int outInSineInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inExponentInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static int outInExponentInterpolation(float startvalue, float endvalue, float time);

    /**
     * interpolates between the startvalue and the endvalue
     * interpolation function like inCircInterpolation, but first easing out until x=0.5 then easing in
     * timeinterval 0-1
     */
    static int outInCircInterpolation(float startvalue, float endvalue, float time);

};

} // namespace voreen
#endif
