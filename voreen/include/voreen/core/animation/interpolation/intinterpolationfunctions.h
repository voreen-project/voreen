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

#ifndef VRN_INTINTERPOLATIONFUNCTIONS_H
#define VRN_INTINTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<int>;
#endif

/**
 * This class offers an interpolation function for int-values. Interpolation: focus on startvalue.
 */
class VRN_CORE_API IntStartInterpolationFunction : public InterpolationFunction<int> {
public:
    IntStartInterpolationFunction();
    virtual std::string getClassName() const { return "IntStartInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: focus on endvalue.
 */
class VRN_CORE_API IntEndInterpolationFunction : public InterpolationFunction<int> {
public:
    IntEndInterpolationFunction();
    virtual std::string getClassName() const { return "IntEndInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: bisection.
 */
class VRN_CORE_API IntStartEndInterpolationFunction : public InterpolationFunction<int> {
public:
    IntStartEndInterpolationFunction();
    virtual std::string getClassName() const { return "IntStartEndInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: linear.
 */
class VRN_CORE_API IntLinearInterpolationFunction : public InterpolationFunction<int> {
public:
    IntLinearInterpolationFunction();
    virtual std::string getClassName() const { return "IntLinearInterpolationFunction"; }
    virtual InterpolationFunction<int>* create() const;
    virtual int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: Catmull-Rom spline.
 */
class VRN_CORE_API IntCatmullRomInterpolationFunction : public MultiPointInterpolationFunction<int> {
public:
    IntCatmullRomInterpolationFunction();
    virtual std::string getClassName() const { return "IntCatmullRomInterpolationFunction"; }
    MultiPointInterpolationFunction<int>* create() const;
    int interpolate(std::vector<PropertyKeyValue<int>*> controlpoints, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quadratic (easing in).
 */
class VRN_CORE_API IntInQuadInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IntInQuadInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: cubicular (easing in).
 */
class VRN_CORE_API IntInCubicInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IntInCubicInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quartetic (easing in).
 */
class VRN_CORE_API IntInQuartInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IntInQuartInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quintic (easing in).
 */
class VRN_CORE_API IntInQuintInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IntInQuintInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: sineousidal (easing in).
 */
class VRN_CORE_API IntInSineInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInSineInterpolationFunction();
    virtual std::string getClassName() const { return "IntInSineInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: exponential (easing in).
 */
class VRN_CORE_API IntInExponentInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IntInExponentInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: circular (easing in).
 */
class VRN_CORE_API IntInCircInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInCircInterpolationFunction();
    virtual std::string getClassName() const { return "IntInCircInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quadratic (easing out).
 */
class VRN_CORE_API IntOutQuadInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutQuadInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: cubicular (easing out).
 */
class VRN_CORE_API IntOutCubicInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutCubicInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quartetic (easing out).
 */
class VRN_CORE_API IntOutQuartInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutQuartInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quintic (easing out).
 */
class VRN_CORE_API IntOutQuintInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutQuintInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: sineousidal (easing out).
 */
class VRN_CORE_API IntOutSineInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutSineInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutSineInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: exponential (easing out).
 */
class VRN_CORE_API IntOutExponentInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutExponentInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: circular (easing out).
 */
class VRN_CORE_API IntOutCircInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutCircInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutCircInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quadratic (easing in, then easing out).
 */
class VRN_CORE_API IntInOutQuadInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInOutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IntInOutQuadInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: cubicular (easing in, then easing out).
 */
class VRN_CORE_API IntInOutCubicInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInOutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IntInOutCubicInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quartetic (easing in, then easing out).
 */
class VRN_CORE_API IntInOutQuartInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInOutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IntInOutQuartInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quintic (easing in, then easing out).
 */
class VRN_CORE_API IntInOutQuintInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInOutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IntInOutQuintInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: sineousidal (easing in, then easing out).
 */
class VRN_CORE_API IntInOutSineInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInOutSineInterpolationFunction();
    virtual std::string getClassName() const { return "IntInOutSineInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: exponential (easing in, then easing out).
 */
class VRN_CORE_API IntInOutExponentInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInOutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IntInOutExponentInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: circular (easing in, then easing out).
 */
class VRN_CORE_API IntInOutCircInterpolationFunction : public InterpolationFunction<int> {
public:
    IntInOutCircInterpolationFunction();
    virtual std::string getClassName() const { return "IntInOutCircInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quadratic (easing out, then easing in).
 */
class VRN_CORE_API IntOutInQuadInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutInQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutInQuadInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: cubicular (easing out, then easing in).
 */
class VRN_CORE_API IntOutInCubicInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutInCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutInCubicInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quartetic (easing out, then easing in).
 */
class VRN_CORE_API IntOutInQuartInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutInQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutInQuartInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: quintic (easing out, then easing in).
 */
class VRN_CORE_API IntOutInQuintInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutInQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutInQuintInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: sineousidal (easing out, then easing in).
 */
class VRN_CORE_API IntOutInSineInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutInSineInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutInSineInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: exponential (easing out, then easing in).
 */
class VRN_CORE_API IntOutInExponentInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutInExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutInExponentInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for int-values. Interpolation: circular (easing out, then easing in).
 */
class VRN_CORE_API IntOutInCircInterpolationFunction : public InterpolationFunction<int> {
public:
    IntOutInCircInterpolationFunction();
    virtual std::string getClassName() const { return "IntOutInCircInterpolationFunction"; }
    InterpolationFunction<int>* create() const;
    int interpolate(int startvalue, int endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

} // namespace voreen

#endif
