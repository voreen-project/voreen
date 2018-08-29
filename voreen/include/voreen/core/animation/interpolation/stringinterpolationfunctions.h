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

#ifndef VRN_VOLUMESTRINGINTERPOLATIONFUNCTIONS_H
#define VRN_VOLUMESTRINGINTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<std::string>;
#endif

/**
 * This class offers an interpolation function for strings. Interpolation: focus on startvalue.
 */
class VRN_CORE_API StringStartInterpolationFunction : public InterpolationFunction<std::string> {
public:
    StringStartInterpolationFunction();
    virtual std::string getClassName() const { return "StringStartInterpolationFunction"; }
    InterpolationFunction<std::string>* create() const;
    std::string interpolate(std::string startvalue, std::string endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for strings. Interpolation: focus on endvalue.
 */
class VRN_CORE_API StringEndInterpolationFunction : public InterpolationFunction<std::string> {
public:
    StringEndInterpolationFunction();
    virtual std::string getClassName() const { return "StringEndInterpolationFunction"; }
    InterpolationFunction<std::string>* create() const;
    std::string interpolate(std::string startvalue, std::string endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for strings. Interpolation: bisection.
 */
class VRN_CORE_API StringStartEndInterpolationFunction : public InterpolationFunction<std::string> {
public:
    StringStartEndInterpolationFunction();
    virtual std::string getClassName() const { return "StringStartEndInterpolationFunction"; }
    InterpolationFunction<std::string>* create() const;
    std::string interpolate(std::string startvalue, std::string endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

} // namespace voreen
#endif
