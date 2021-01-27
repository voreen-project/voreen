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

#ifndef VRN_MAT2INTERPOLATIONFUNCTIONS_H
#define VRN_MAT2INTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<tgt::mat2>;
#endif

/**
 * This class offers an interpolation function for mat2-values. Interpolation: focus on startvalue.
 */
class VRN_CORE_API Mat2StartInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2StartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2StartInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2-values. Interpolation: focus on endvalue.
 */
class VRN_CORE_API Mat2EndInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2EndInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2EndInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2-values. Interpolation: bisection.
 */
class VRN_CORE_API Mat2StartEndInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2StartEndInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2StartEndInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: linear.
 */
class VRN_CORE_API Mat2LinearInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2LinearInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2LinearInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quadratic (easing in).
 */
class VRN_CORE_API Mat2InQuadInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: cubicular (easing in).
 */
class VRN_CORE_API Mat2InCubicInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quartetic (easing in).
 */
class VRN_CORE_API Mat2InQuartInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quintic (easing in).
 */
class VRN_CORE_API Mat2InQuintInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: sineousidal (easing in).
 */
class VRN_CORE_API Mat2InSineInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: exponential (easing in).
 */
class VRN_CORE_API Mat2InExponentInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: circular (easing in).
 */
class VRN_CORE_API Mat2InCircInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quadratic (easing out).
 */
class VRN_CORE_API Mat2OutQuadInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: cubicular (easing out).
 */
class VRN_CORE_API Mat2OutCubicInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quartetic (easing out).
 */
class VRN_CORE_API Mat2OutQuartInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quintic (easing out).
 */
class VRN_CORE_API Mat2OutQuintInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: sineousidal (easing out).
 */
class VRN_CORE_API Mat2OutSineInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: exponential (easing out).
 */
class VRN_CORE_API Mat2OutExponentInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: circular (easing out).
 */
class VRN_CORE_API Mat2OutCircInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quadratic (easing in, then easing out).
 */
class VRN_CORE_API Mat2InOutQuadInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InOutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InOutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: cubicular (easing in, then easing out).
 */
class VRN_CORE_API Mat2InOutCubicInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InOutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InOutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quartetic (easing in, then easing out).
 */
class VRN_CORE_API Mat2InOutQuartInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InOutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InOutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quintic (easing in, then easing out).
 */
class VRN_CORE_API Mat2InOutQuintInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InOutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InOutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: sineousidal (easing in, then easing out).
 */
class VRN_CORE_API Mat2InOutSineInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InOutSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InOutSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: exponential (easing in, then easing out).
 */
class VRN_CORE_API Mat2InOutExponentInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InOutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InOutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: circular (easing in, then easing out).
 */
class VRN_CORE_API Mat2InOutCircInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2InOutCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2InOutCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quadratic (easing out, then easing in).
 */
class VRN_CORE_API Mat2OutInQuadInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutInQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutInQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: cubicular (easing out, then easing in).
 */
class VRN_CORE_API Mat2OutInCubicInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutInCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutInCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quartetic (easing out, then easing in).
 */
class VRN_CORE_API Mat2OutInQuartInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutInQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutInQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: quintic (easing out, then easing in).
 */
class VRN_CORE_API Mat2OutInQuintInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutInQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutInQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: sineousidal (easing out, then easing in).
 */
class VRN_CORE_API Mat2OutInSineInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutInSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutInSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: exponential (easing out, then easing in).
 */
class VRN_CORE_API Mat2OutInExponentInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutInExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutInExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat2. Interpolation: circular (easing out, then easing in).
 */
class VRN_CORE_API Mat2OutInCircInterpolationFunction : public InterpolationFunction<tgt::mat2> {
public:
    Mat2OutInCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat2OutInCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat2>* create() const;
    tgt::mat2 interpolate(tgt::mat2 startvalue, tgt::mat2 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

} // namespace voreen
#endif
