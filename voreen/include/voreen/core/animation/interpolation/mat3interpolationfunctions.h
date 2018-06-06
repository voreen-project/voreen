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

#ifndef VRN_MAT3INTERPOLATIONFUNCTIONS_H
#define VRN_MAT3INTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<tgt::mat3>;
#endif

/**
 * This class offers an interpolation function for mat3-values. Interpolation: focus on startvalue.
 */
class VRN_CORE_API Mat3StartInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3StartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3StartInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3-values. Interpolation: focus on endvalue.
 */
class VRN_CORE_API Mat3EndInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3EndInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3EndInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3-values. Interpolation: bisection.
 */
class VRN_CORE_API Mat3StartEndInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3StartEndInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3StartEndInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: linear.
 */
class VRN_CORE_API Mat3LinearInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3LinearInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3LinearInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quadratic (easing in).
 */
class VRN_CORE_API Mat3InQuadInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: cubicular (easing in).
 */
class VRN_CORE_API Mat3InCubicInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quartetic (easing in).
 */
class VRN_CORE_API Mat3InQuartInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quintic (easing in).
 */
class VRN_CORE_API Mat3InQuintInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: sineousidal (easing in).
 */
class VRN_CORE_API Mat3InSineInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: exponential (easing in).
 */
class VRN_CORE_API Mat3InExponentInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: circular (easing in).
 */
class VRN_CORE_API Mat3InCircInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quadratic (easing out).
 */
class VRN_CORE_API Mat3OutQuadInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: cubicular (easing out).
 */
class VRN_CORE_API Mat3OutCubicInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quartetic (easing out).
 */
class VRN_CORE_API Mat3OutQuartInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quintic (easing out).
 */
class VRN_CORE_API Mat3OutQuintInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: sineousidal (easing out).
 */
class VRN_CORE_API Mat3OutSineInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: exponential (easing out).
 */
class VRN_CORE_API Mat3OutExponentInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: circular (easing out).
 */
class VRN_CORE_API Mat3OutCircInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quadratic (easing in, then easing out).
 */
class VRN_CORE_API Mat3InOutQuadInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InOutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InOutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: cubicular (easing in, then easing out).
 */
class VRN_CORE_API Mat3InOutCubicInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InOutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InOutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quartetic (easing in, then easing out).
 */
class VRN_CORE_API Mat3InOutQuartInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InOutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InOutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quintic (easing in, then easing out).
 */
class VRN_CORE_API Mat3InOutQuintInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InOutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InOutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: sineousidal (easing in, then easing out).
 */
class VRN_CORE_API Mat3InOutSineInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InOutSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InOutSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: exponential (easing in, then easing out).
 */
class VRN_CORE_API Mat3InOutExponentInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InOutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InOutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: circular (easing in, then easing out).
 */
class VRN_CORE_API Mat3InOutCircInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3InOutCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3InOutCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quadratic (easing out, then easing in).
 */
class VRN_CORE_API Mat3OutInQuadInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutInQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutInQuadInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: cubicular (easing out, then easing in).
 */
class VRN_CORE_API Mat3OutInCubicInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutInCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutInCubicInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quartetic (easing out, then easing in).
 */
class VRN_CORE_API Mat3OutInQuartInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutInQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutInQuartInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: quintic (easing out, then easing in).
 */
class VRN_CORE_API Mat3OutInQuintInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutInQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutInQuintInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: sineousidal (easing out, then easing in).
 */
class VRN_CORE_API Mat3OutInSineInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutInSineInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutInSineInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: exponential (easing out, then easing in).
 */
class VRN_CORE_API Mat3OutInExponentInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutInExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutInExponentInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for mat3. Interpolation: circular (easing out, then easing in).
 */
class VRN_CORE_API Mat3OutInCircInterpolationFunction : public InterpolationFunction<tgt::mat3> {
public:
    Mat3OutInCircInterpolationFunction();
    virtual std::string getClassName() const { return "Mat3OutInCircInterpolationFunction"; }
    InterpolationFunction<tgt::mat3>* create() const;
    tgt::mat3 interpolate(tgt::mat3 startvalue, tgt::mat3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

} // namespace voreen
#endif
