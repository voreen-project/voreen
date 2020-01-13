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

#ifndef VRN_MAT4INTERPOLATIONFUNCTIONS_H
#define VRN_MAT4INTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<tgt::mat4>;
#endif

/**
 * This class offers an interpolation function for mat4-values. Interpolation: focus on startvalue.
 */
class VRN_CORE_API Mat4StartInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4StartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4StartInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4-values. Interpolation: focus on endvalue.
 */
class VRN_CORE_API Mat4EndInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4EndInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4EndInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4-values. Interpolation: bisection.
 */
class VRN_CORE_API Mat4StartEndInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4StartEndInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4StartEndInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: linear.
 */
class VRN_CORE_API Mat4LinearInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4LinearInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4LinearInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quadratic (easing in).
 */
class VRN_CORE_API Mat4InQuadInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: cubicular (easing in).
 */
class VRN_CORE_API Mat4InCubicInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quartetic (easing in).
 */
class VRN_CORE_API Mat4InQuartInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quintic (easing in).
 */
class VRN_CORE_API Mat4InQuintInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: sineousidal (easing in).
 */
class VRN_CORE_API Mat4InSineInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: exponential (easing in).
 */
class VRN_CORE_API Mat4InExponentInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: circular (easing in).
 */
class VRN_CORE_API Mat4InCircInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quadratic (easing out).
 */
class VRN_CORE_API Mat4OutQuadInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: cubicular (easing out).
 */
class VRN_CORE_API Mat4OutCubicInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quartetic (easing out).
 */
class VRN_CORE_API Mat4OutQuartInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quintic (easing out).
 */
class VRN_CORE_API Mat4OutQuintInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: sineousidal (easing out).
 */
class VRN_CORE_API Mat4OutSineInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: exponential (easing out).
 */
class VRN_CORE_API Mat4OutExponentInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: circular (easing out).
 */
class VRN_CORE_API Mat4OutCircInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quadratic (easing in, then easing out).
 */
class VRN_CORE_API Mat4InOutQuadInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InOutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InOutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: cubicular (easing in, then easing out).
 */
class VRN_CORE_API Mat4InOutCubicInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InOutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InOutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quartetic (easing in, then easing out).
 */
class VRN_CORE_API Mat4InOutQuartInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InOutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InOutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quintic (easing in, then easing out).
 */
class VRN_CORE_API Mat4InOutQuintInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InOutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InOutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: sineousidal (easing in, then easing out).
 */
class VRN_CORE_API Mat4InOutSineInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InOutSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InOutSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: exponential (easing in, then easing out).
 */
class VRN_CORE_API Mat4InOutExponentInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InOutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InOutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: circular (easing in, then easing out).
 */
class VRN_CORE_API Mat4InOutCircInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4InOutCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4InOutCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quadratic (easing out, then easing in).
 */
class VRN_CORE_API Mat4OutInQuadInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutInQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutInQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: cubicular (easing out, then easing in).
 */
class VRN_CORE_API Mat4OutInCubicInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutInCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutInCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quartetic (easing out, then easing in).
 */
class VRN_CORE_API Mat4OutInQuartInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutInQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutInQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: quintic (easing out, then easing in).
 */
class VRN_CORE_API Mat4OutInQuintInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutInQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutInQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: sineousidal (easing out, then easing in).
 */
class VRN_CORE_API Mat4OutInSineInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutInSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutInSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: exponential (easing out, then easing in).
 */
class VRN_CORE_API Mat4OutInExponentInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutInExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutInExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat4. Interpolation: circular (easing out, then easing in).
 */
class VRN_CORE_API Mat4OutInCircInterpolationFunction : public InterpolationFunction<tgt::mat4> {
public:
    Mat4OutInCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat4OutInCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat4>* create() const;
    tgt::mat4 interpolate(tgt::mat4 startvalue, tgt::mat4 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

} // namespace voreen
#endif
