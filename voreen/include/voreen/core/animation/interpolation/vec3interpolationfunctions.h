/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_VEC3INTERPOLATIONFUNCTIONS_H
#define VRN_VEC3INTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<tgt::vec3>;
#endif

/**
 * This class offers an interpolation function for vec3-values. Interpolation: focus on startvalue.
 */
class VRN_CORE_API Vec3StartInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3StartInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3StartInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3-values. Interpolation: focus on endvalue.
 */
class VRN_CORE_API Vec3EndInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3EndInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3EndInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3-values. Interpolation: bisection.
 */
class VRN_CORE_API Vec3StartEndInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3StartEndInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3StartEndInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: linear.
 */
class VRN_CORE_API Vec3LinearInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3LinearInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3LinearInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: spherical linear.
 */
class VRN_CORE_API Vec3SphericalLinearInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3SphericalLinearInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3SphericalLinearInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: spherical cubic.
 */
class VRN_CORE_API Vec3SphericalCubicInterpolationFunction : public MultiPointInterpolationFunction<tgt::vec3> {
public:
    Vec3SphericalCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3SphericalCubicInterpolationFunction"; }
    MultiPointInterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(std::vector<PropertyKeyValue<tgt::vec3>*> controlpoints, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: Catmull-Rom spline.
 */
class VRN_CORE_API Vec3CatmullRomInterpolationFunction : public MultiPointInterpolationFunction<tgt::vec3>{
public:
    Vec3CatmullRomInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3CatmullRomInterpolationFunction"; }
    MultiPointInterpolationFunction<tgt::vec3>* create() const;

    tgt::vec3 interpolate(std::vector<PropertyKeyValue<tgt::vec3>*> controlpoints, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};


/**
 * This class offers an interpolation function for vec3. Interpolation: quadratic (easing in).
 */
class VRN_CORE_API Vec3InQuadInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InQuadInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: cubicular (easing in).
 */
class VRN_CORE_API Vec3InCubicInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InCubicInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: quartetic (easing in).
 */
class VRN_CORE_API Vec3InQuartInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InQuartInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: quintic (easing in).
 */
class VRN_CORE_API Vec3InQuintInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InQuintInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: sineousidal (easing in).
 */
class VRN_CORE_API Vec3InSineInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InSineInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InSineInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: exponential (easing in).
 */
class VRN_CORE_API Vec3InExponentInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InExponentInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: circular (easing in).
 */
class VRN_CORE_API Vec3InCircInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InCircInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InCircInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: quadratic (easing out).
 */
class VRN_CORE_API Vec3OutQuadInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: cubicular (easing out).
 */
class VRN_CORE_API Vec3OutCubicInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: quartetic (easing out).
 */
class VRN_CORE_API Vec3OutQuartInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: quintic (easing out).
 */
class VRN_CORE_API Vec3OutQuintInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: sineousidal (easing out).
 */
class VRN_CORE_API Vec3OutSineInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutSineInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutSineInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: exponential (easing out).
 */
class VRN_CORE_API Vec3OutExponentInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: circular (easing out).
 */
class VRN_CORE_API Vec3OutCircInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutCircInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutCircInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: quadratic (easing in, then easing out).
 */
class VRN_CORE_API Vec3InOutQuadInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InOutQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InOutQuadInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: cubicular (easing in, then easing out).
 */
class VRN_CORE_API Vec3InOutCubicInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InOutCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InOutCubicInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: quartetic (easing in, then easing out).
 */
class VRN_CORE_API Vec3InOutQuartInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InOutQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InOutQuartInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: quintic (easing in, then easing out).
 */
class VRN_CORE_API Vec3InOutQuintInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InOutQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InOutQuintInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: sineousidal (easing in, then easing out).
 */
class VRN_CORE_API Vec3InOutSineInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InOutSineInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InOutSineInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: exponential (easing in, then easing out).
 */
class VRN_CORE_API Vec3InOutExponentInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InOutExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InOutExponentInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: circular (easing in, then easing out).
 */
class VRN_CORE_API Vec3InOutCircInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3InOutCircInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3InOutCircInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: quadratic (easing out, then easing in).
 */
class VRN_CORE_API Vec3OutInQuadInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutInQuadInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutInQuadInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: cubicular (easing out, then easing in).
 */
class VRN_CORE_API Vec3OutInCubicInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutInCubicInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutInCubicInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: quartetic (easing out, then easing in).
 */
class VRN_CORE_API Vec3OutInQuartInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutInQuartInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutInQuartInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: quintic (easing out, then easing in).
 */
class VRN_CORE_API Vec3OutInQuintInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutInQuintInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutInQuintInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: sineousidal (easing out, then easing in).
 */
class VRN_CORE_API Vec3OutInSineInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutInSineInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutInSineInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: exponential (easing out, then easing in).
 */
class VRN_CORE_API Vec3OutInExponentInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutInExponentInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutInExponentInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for vec3. Interpolation: circular (easing out, then easing in).
 */
class VRN_CORE_API Vec3OutInCircInterpolationFunction : public InterpolationFunction<tgt::vec3> {
public:
    Vec3OutInCircInterpolationFunction();
    virtual std::string getClassName() const { return "Vec3OutInCircInterpolationFunction"; }
    InterpolationFunction<tgt::vec3>* create() const;
    tgt::vec3 interpolate(tgt::vec3 startvalue, tgt::vec3 endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

} // namespace voreen
#endif
