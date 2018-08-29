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

#ifndef VRN_IVEC3INTERPOLATIONFUNCTIONS_H
#define VRN_IVEC3INTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<tgt::ivec3>;
#endif

/**
 * This class offers an interpolation function for ivec3-values. Interpolation: focus on startvalue.
 */
class VRN_CORE_API IVec3StartInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3StartInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3StartInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3-values. Interpolation: focus on endvalue.
 */
class VRN_CORE_API IVec3EndInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3EndInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3EndInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3-values. Interpolation: bisection.
 */
class VRN_CORE_API IVec3StartEndInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3StartEndInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3StartEndInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: linear.
 */
class VRN_CORE_API IVec3LinearInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3LinearInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3LinearInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quadratic (easing in).
 */
class VRN_CORE_API IVec3InQuadInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InQuadInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: cubicular (easing in).
 */
class VRN_CORE_API IVec3InCubicInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InCubicInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quartetic (easing in).
 */
class VRN_CORE_API IVec3InQuartInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InQuartInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quintic (easing in).
 */
class VRN_CORE_API IVec3InQuintInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InQuintInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: sineousidal (easing in).
 */
class VRN_CORE_API IVec3InSineInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InSineInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InSineInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: exponential (easing in).
 */
class VRN_CORE_API IVec3InExponentInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InExponentInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: circular (easing in).
 */
class VRN_CORE_API IVec3InCircInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InCircInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InCircInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quadratic (easing out).
 */
class VRN_CORE_API IVec3OutQuadInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: cubicular (easing out).
 */
class VRN_CORE_API IVec3OutCubicInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quartetic (easing out).
 */
class VRN_CORE_API IVec3OutQuartInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quintic (easing out).
 */
class VRN_CORE_API IVec3OutQuintInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: sineousidal (easing out).
 */
class VRN_CORE_API IVec3OutSineInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutSineInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutSineInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: exponential (easing out).
 */
class VRN_CORE_API IVec3OutExponentInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: circular (easing out).
 */
class VRN_CORE_API IVec3OutCircInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutCircInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutCircInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quadratic (easing in, then easing out).
 */
class VRN_CORE_API IVec3InOutQuadInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InOutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InOutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: cubicular (easing in, then easing out).
 */
class VRN_CORE_API IVec3InOutCubicInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InOutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InOutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quartetic (easing in, then easing out).
 */
class VRN_CORE_API IVec3InOutQuartInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InOutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InOutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quintic (easing in, then easing out).
 */
class VRN_CORE_API IVec3InOutQuintInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InOutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InOutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: sineousidal (easing in, then easing out).
 */
class VRN_CORE_API IVec3InOutSineInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InOutSineInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InOutSineInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: exponential (easing in, then easing out).
 */
class VRN_CORE_API IVec3InOutExponentInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InOutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InOutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: circular (easing in, then easing out).
 */
class VRN_CORE_API IVec3InOutCircInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3InOutCircInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3InOutCircInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quadratic (easing out, then easing in).
 */
class VRN_CORE_API IVec3OutInQuadInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutInQuadInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutInQuadInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: cubicular (easing out, then easing in).
 */
class VRN_CORE_API IVec3OutInCubicInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutInCubicInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutInCubicInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quartetic (easing out, then easing in).
 */
class VRN_CORE_API IVec3OutInQuartInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutInQuartInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutInQuartInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: quintic (easing out, then easing in).
 */
class VRN_CORE_API IVec3OutInQuintInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutInQuintInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutInQuintInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: sineousidal (easing out, then easing in).
 */
class VRN_CORE_API IVec3OutInSineInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutInSineInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutInSineInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: exponential (easing out, then easing in).
 */
class VRN_CORE_API IVec3OutInExponentInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutInExponentInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutInExponentInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for ivec3. Interpolation: circular (easing out, then easing in).
 */
class VRN_CORE_API IVec3OutInCircInterpolationFunction : public InterpolationFunction<tgt::ivec3> {
public:
    IVec3OutInCircInterpolationFunction();
    virtual std::string getClassName() const { return "IVec3OutInCircInterpolationFunction"; }
    InterpolationFunction<tgt::ivec3>* create() const;
    tgt::ivec3 interpolate(tgt::ivec3 startvalue, tgt::ivec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

} // namespace voreen
#endif
