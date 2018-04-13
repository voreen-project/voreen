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

#ifndef VRN_VEC2INTERPOLATIONFUNCTIONS_H
#define VRN_VEC2INTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<tgt::vec2>;
#endif

/**
 * This class offers an interpolation function for vec2-values. Interpolation: focus on startvalue.
 */
class VRN_CORE_API Vec2StartInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2StartInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2StartInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2-values. Interpolation: focus on endvalue.
 */
class VRN_CORE_API Vec2EndInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2EndInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2EndInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2-values. Interpolation: bisection.
 */
class VRN_CORE_API Vec2StartEndInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2StartEndInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2StartEndInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: linear.
 */
class VRN_CORE_API Vec2LinearInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2LinearInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2LinearInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: spherical linear.
 */
class VRN_CORE_API Vec2SphericalLinearInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2SphericalLinearInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2SphericalLinearInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quadratic (easing in).
 */
class VRN_CORE_API Vec2InQuadInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InQuadInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: cubicular (easing in).
 */
class VRN_CORE_API Vec2InCubicInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InCubicInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quartetic (easing in).
 */
class VRN_CORE_API Vec2InQuartInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InQuartInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quintic (easing in).
 */
class VRN_CORE_API Vec2InQuintInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InQuintInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: sineousidal (easing in).
 */
class VRN_CORE_API Vec2InSineInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InSineInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InSineInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: exponential (easing in).
 */
class VRN_CORE_API Vec2InExponentInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InExponentInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: circular (easing in).
 */
class VRN_CORE_API Vec2InCircInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InCircInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InCircInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quadratic (easing out).
 */
class VRN_CORE_API Vec2OutQuadInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: cubicular (easing out).
 */
class VRN_CORE_API Vec2OutCubicInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quartetic (easing out).
 */
class VRN_CORE_API Vec2OutQuartInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quintic (easing out).
 */
class VRN_CORE_API Vec2OutQuintInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: sineousidal (easing out).
 */
class VRN_CORE_API Vec2OutSineInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutSineInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutSineInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: exponential (easing out).
 */
class VRN_CORE_API Vec2OutExponentInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: circular (easing out).
 */
class VRN_CORE_API Vec2OutCircInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutCircInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutCircInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quadratic (easing in, then easing out).
 */
class VRN_CORE_API Vec2InOutQuadInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InOutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InOutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: cubicular (easing in, then easing out).
 */
class VRN_CORE_API Vec2InOutCubicInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InOutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InOutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quartetic (easing in, then easing out).
 */
class VRN_CORE_API Vec2InOutQuartInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InOutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InOutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quintic (easing in, then easing out).
 */
class VRN_CORE_API Vec2InOutQuintInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InOutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InOutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: sineousidal (easing in, then easing out).
 */
class VRN_CORE_API Vec2InOutSineInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InOutSineInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InOutSineInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: exponential (easing in, then easing out).
 */
class VRN_CORE_API Vec2InOutExponentInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InOutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InOutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: circular (easing in, then easing out).
 */
class VRN_CORE_API Vec2InOutCircInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2InOutCircInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2InOutCircInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quadratic (easing out, then easing in).
 */
class VRN_CORE_API Vec2OutInQuadInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutInQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutInQuadInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: cubicular (easing out, then easing in).
 */
class VRN_CORE_API Vec2OutInCubicInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutInCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutInCubicInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quartetic (easing out, then easing in).
 */
class VRN_CORE_API Vec2OutInQuartInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutInQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutInQuartInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: quintic (easing out, then easing in).
 */
class VRN_CORE_API Vec2OutInQuintInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutInQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutInQuintInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: sineousidal (easing out, then easing in).
 */
class VRN_CORE_API Vec2OutInSineInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutInSineInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutInSineInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: exponential (easing out, then easing in).
 */
class VRN_CORE_API Vec2OutInExponentInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutInExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutInExponentInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec2. Interpolation: circular (easing out, then easing in).
 */
class VRN_CORE_API Vec2OutInCircInterpolationFunction : public InterpolationFunction<tgt::vec2> {
public:
    Vec2OutInCircInterpolationFunction();
    virtual std::string getClassName() const { return "Vec2OutInCircInterpolationFunction"; }
    InterpolationFunction<tgt::vec2>* create() const;
    tgt::vec2 interpolate(tgt::vec2 startvalue, tgt::vec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

} // namespace voreen

#endif
