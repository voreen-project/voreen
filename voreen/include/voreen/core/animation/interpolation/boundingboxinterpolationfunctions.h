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

#ifndef VRN_BOUNDINGBOXINTERPOLATIONFUNCTIONS_H
#define VRN_BOUNDINGBOXINTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"
#include "tgt/bounds.h"


namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<tgt::IntBounds>;
#endif

/**
 * This class offers an interpolation function for boundingbox-values. Interpolation: linear.
 */
class VRN_CORE_API IntBoundingboxLinearInterpolationFunction : public InterpolationFunction<tgt::IntBounds> {
public:
    IntBoundingboxLinearInterpolationFunction(){}
    virtual std::string getClassName() const { return "IntBoundingboxLinearInterpolationFunction"; }
    InterpolationFunction<tgt::IntBounds>* create() const {return new IntBoundingboxLinearInterpolationFunction;}
    tgt::IntBounds interpolate(tgt::IntBounds startvalue, tgt::IntBounds endvalue, float time) const;
    std::string getGuiName() const { return "linear interpolation";}
    std::string getCategory() const{ return "boundingbox"; }
};


/**
 * This class offers an interpolation function for boundingbox-values. Interpolation: focus on startvalue.
 */
class VRN_CORE_API IntBoundingboxStartInterpolationFunction : public InterpolationFunction<tgt::IntBounds> {
public:
    IntBoundingboxStartInterpolationFunction(){}
    virtual std::string getClassName() const { return "IntBoundingboxStartInterpolationFunction"; }
    InterpolationFunction<tgt::IntBounds>* create() const {return new IntBoundingboxStartInterpolationFunction;}
    tgt::IntBounds interpolate(tgt::IntBounds startvalue, tgt::IntBounds endvalue, float time) const;
    std::string getGuiName() const { return "focus on startvalue";}
    std::string getCategory() const{ return "boundingbox"; }
};


/**
 * This class offers an interpolation function for boundingbox-values. Interpolation: focus on endvalue.
 */
class VRN_CORE_API IntBoundingboxEndInterpolationFunction : public InterpolationFunction<tgt::IntBounds> {
public:
    IntBoundingboxEndInterpolationFunction(){}
    virtual std::string getClassName() const { return "IntBoundingboxEndInterpolationFunction"; }
    InterpolationFunction<tgt::IntBounds>* create() const {return new IntBoundingboxEndInterpolationFunction;}
    tgt::IntBounds interpolate(tgt::IntBounds startvalue, tgt::IntBounds endvalue, float time) const;
    std::string getGuiName() const { return "focus on endvalue";}
    std::string getCategory() const{ return "boundingbox"; }
};

/**
 * This class offers an interpolation function for boundingbox-values. Interpolation: bisection.
 */
class VRN_CORE_API IntBoundingboxStartEndInterpolationFunction : public InterpolationFunction<tgt::IntBounds> {
public:
    IntBoundingboxStartEndInterpolationFunction(){}
    virtual std::string getClassName() const { return "IntBoundingboxStartEndInterpolationFunction"; }
    InterpolationFunction<tgt::IntBounds>* create() const {return new IntBoundingboxStartEndInterpolationFunction;}
    tgt::IntBounds interpolate(tgt::IntBounds startvalue, tgt::IntBounds endvalue, float time) const;
    std::string getGuiName() const { return "bisection";}
    std::string getCategory() const{ return "boundingbox"; }
};


/**
 * This class offers an interpolation function for boundingbox-values. Interpolation: linear.
 */
class VRN_CORE_API FloatBoundingboxLinearInterpolationFunction : public InterpolationFunction<tgt::Bounds> {
public:
    FloatBoundingboxLinearInterpolationFunction(){}
    virtual std::string getClassName() const { return "FloatBoundingboxLinearInterpolationFunction"; }
    InterpolationFunction<tgt::Bounds>* create() const {return new FloatBoundingboxLinearInterpolationFunction;}
    tgt::Bounds interpolate(tgt::Bounds startvalue, tgt::Bounds endvalue, float time) const;
    std::string getGuiName() const { return "linear interpolation";}
    std::string getCategory() const{ return "boundingbox"; }
};


/**
 * This class offers an interpolation function for boundingbox-values. Interpolation: focus on startvalue.
 */
class VRN_CORE_API FloatBoundingboxStartInterpolationFunction : public InterpolationFunction<tgt::Bounds> {
public:
    FloatBoundingboxStartInterpolationFunction(){}
    virtual std::string getClassName() const { return "FloatBoundingboxStartInterpolationFunction"; }
    InterpolationFunction<tgt::Bounds>* create() const {return new FloatBoundingboxStartInterpolationFunction;}
    tgt::Bounds interpolate(tgt::Bounds startvalue, tgt::Bounds endvalue, float time) const;
    std::string getGuiName() const { return "focus on startvalue";}
    std::string getCategory() const{ return "boundingbox"; }
};


/**
 * This class offers an interpolation function for boundingbox-values. Interpolation: focus on endvalue.
 */
class VRN_CORE_API FloatBoundingboxEndInterpolationFunction : public InterpolationFunction<tgt::Bounds> {
public:
    FloatBoundingboxEndInterpolationFunction(){}
    virtual std::string getClassName() const { return "FloatBoundingboxEndInterpolationFunction"; }
    InterpolationFunction<tgt::Bounds>* create() const {return new FloatBoundingboxEndInterpolationFunction;}
    tgt::Bounds interpolate(tgt::Bounds startvalue, tgt::Bounds endvalue, float time) const;
    std::string getGuiName() const { return "focus on endvalue";}
    std::string getCategory() const{ return "boundingbox"; }
};

/**
 * This class offers an interpolation function for boundingbox-values. Interpolation: bisection.
 */
class VRN_CORE_API FloatBoundingboxStartEndInterpolationFunction : public InterpolationFunction<tgt::Bounds> {
public:
    FloatBoundingboxStartEndInterpolationFunction(){}
    virtual std::string getClassName() const { return "FloatBoundingboxStartEndInterpolationFunction"; }
    InterpolationFunction<tgt::Bounds>* create() const {return new FloatBoundingboxStartEndInterpolationFunction;}
    tgt::Bounds interpolate(tgt::Bounds startvalue, tgt::Bounds endvalue, float time) const;
    std::string getGuiName() const { return "bisection";}
    std::string getCategory() const{ return "boundingbox"; }
};

} // namespace voreen

#endif
