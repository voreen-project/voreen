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

#ifndef VRN_IVEC2INTERPOLATIONFUNCTIONS_H
#define VRN_IVEC2INTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"
#include "tgt/vector.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<tgt::ivec2>;
#endif

/**
 * This class offers an interpolation function for ivec2-values. Interpolation: focus on startvalue.
 */
class VRN_CORE_API IVec2StartInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2StartInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2StartInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2-values. Interpolation: focus on endvalue.
 */
class VRN_CORE_API IVec2EndInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2EndInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2EndInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2-values. Interpolation: bisection.
 */
class VRN_CORE_API IVec2StartEndInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2StartEndInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2StartEndInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: linear.
 */
class VRN_CORE_API IVec2LinearInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2LinearInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2LinearInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quadratic (easing in).
 */
class VRN_CORE_API IVec2InQuadInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InQuadInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: cubicular (easing in).
 */
class VRN_CORE_API IVec2InCubicInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InCubicInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quartetic (easing in).
 */
class VRN_CORE_API IVec2InQuartInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InQuartInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quintic (easing in).
 */
class VRN_CORE_API IVec2InQuintInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InQuintInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: sineousidal (easing in).
 */
class VRN_CORE_API IVec2InSineInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InSineInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InSineInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: exponential (easing in).
 */
class VRN_CORE_API IVec2InExponentInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InExponentInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: circular (easing in).
 */
class VRN_CORE_API IVec2InCircInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InCircInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InCircInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quadratic (easing out).
 */
class VRN_CORE_API IVec2OutQuadInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: cubicular (easing out).
 */
class VRN_CORE_API IVec2OutCubicInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quartetic (easing out).
 */
class VRN_CORE_API IVec2OutQuartInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quintic (easing out).
 */
class VRN_CORE_API IVec2OutQuintInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: sineousidal (easing out).
 */
class VRN_CORE_API IVec2OutSineInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutSineInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutSineInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: exponential (easing out).
 */
class VRN_CORE_API IVec2OutExponentInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: circular (easing out).
 */
class VRN_CORE_API IVec2OutCircInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutCircInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutCircInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quadratic (easing in, then easing out).
 */
class VRN_CORE_API IVec2InOutQuadInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InOutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InOutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: cubicular (easing in, then easing out).
 */
class VRN_CORE_API IVec2InOutCubicInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InOutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InOutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quartetic (easing in, then easing out).
 */
class VRN_CORE_API IVec2InOutQuartInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InOutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InOutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quintic (easing in, then easing out).
 */
class VRN_CORE_API IVec2InOutQuintInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InOutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InOutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: sineousidal (easing in, then easing out).
 */
class VRN_CORE_API IVec2InOutSineInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InOutSineInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InOutSineInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: exponential (easing in, then easing out).
 */
class VRN_CORE_API IVec2InOutExponentInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InOutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InOutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: circular (easing in, then easing out).
 */
class VRN_CORE_API IVec2InOutCircInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2InOutCircInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2InOutCircInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quadratic (easing out, then easing in).
 */
class VRN_CORE_API IVec2OutInQuadInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutInQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutInQuadInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: cubicular (easing out, then easing in).
 */
class VRN_CORE_API IVec2OutInCubicInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutInCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutInCubicInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quartetic (easing out, then easing in).
 */
class VRN_CORE_API IVec2OutInQuartInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutInQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutInQuartInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: quintic (easing out, then easing in).
 */
class VRN_CORE_API IVec2OutInQuintInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutInQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutInQuintInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: sineousidal (easing out, then easing in).
 */
class VRN_CORE_API IVec2OutInSineInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutInSineInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutInSineInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: exponential (easing out, then easing in).
 */
class VRN_CORE_API IVec2OutInExponentInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutInExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutInExponentInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec2. Interpolation: circular (easing out, then easing in).
 */
class VRN_CORE_API IVec2OutInCircInterpolationFunction : public InterpolationFunction<tgt::ivec2> {
public:
    IVec2OutInCircInterpolationFunction();
    virtual std::string getClassName() const { return "IVec2OutInCircInterpolationFunction"; }
    InterpolationFunction<tgt::ivec2>* create() const;
    tgt::ivec2 interpolate(tgt::ivec2 startvalue, tgt::ivec2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

} // namespace voreen
#endif
