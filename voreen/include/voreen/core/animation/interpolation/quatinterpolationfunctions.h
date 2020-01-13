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

#ifndef VRN_QUATINTERPOLATIONFUNCTIONS_H
#define VRN_QUATINTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"
#include "tgt/quaternion.h"

namespace voreen {

/**
 * This class offers an interpolation function for quat-values. Interpolation: focus on startvalue.
 */
class QuatStartInterpolationFunction : public InterpolationFunction<tgt::quat> {
public:
    QuatStartInterpolationFunction();
    virtual std::string getClassName() const { return "QuatStartInterpolationFunction"; }
    InterpolationFunction<tgt::quat>* create() const;
    tgt::quat interpolate(tgt::quat startvalue, tgt::quat endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for quat-values. Interpolation: focus on endvalue.
 */
class QuatEndInterpolationFunction : public InterpolationFunction<tgt::quat> {
public:
    QuatEndInterpolationFunction();
    virtual std::string getClassName() const { return "QuatEndInterpolationFunction"; }
    InterpolationFunction<tgt::quat>* create() const;
    tgt::quat interpolate(tgt::quat startvalue, tgt::quat endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for quat-values. Interpolation: bisection.
 */
class QuatStartEndInterpolationFunction : public InterpolationFunction<tgt::quat> {
public:
    QuatStartEndInterpolationFunction();
    virtual std::string getClassName() const { return "QuatStartEndInterpolationFunction"; }
    InterpolationFunction<tgt::quat>* create() const;
    tgt::quat interpolate(tgt::quat startvalue, tgt::quat endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for quat. Interpolation: linear.
 */
class QuatLinearInterpolationFunction : public InterpolationFunction<tgt::quat> {
public:
    QuatLinearInterpolationFunction();
    virtual std::string getClassName() const { return "QuatLinearInterpolationFunction"; }
    InterpolationFunction<tgt::quat>* create() const;
    tgt::quat interpolate(tgt::quat startvalue, tgt::quat endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for quat. Interpolation: spherical linear.
 */
class QuatSphericalLinearInterpolationFunction : public InterpolationFunction<tgt::quat> {
public:
    QuatSphericalLinearInterpolationFunction();
    virtual std::string getClassName() const { return "QuatSphericalLinearInterpolationFunction"; }
    InterpolationFunction<tgt::quat>* create() const;
    tgt::quat interpolate(tgt::quat startvalue, tgt::quat endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

} //namespace voreen

#endif
