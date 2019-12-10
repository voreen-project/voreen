/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_TRANSFUNC2DPRIMITIVESINTERPOLATIONFUNCTIONS_H
#define VRN_TRANSFUNC2DPRIMITIVESINTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/transfunc2dprimitives.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<TransFunc2DPrimitives*>;
#endif

/// Base class for transfer function interpolation, provides commonly used methods.
class VRN_CORE_API TransFunc2DPrimitivesInterpolationFunctionBase : public InterpolationFunction<TransFunc2DPrimitives*> {
protected:
    static GLubyte* convertTextureToRGBA(tgt::ivec3 dim, GLubyte* textur, GLuint inputformat);
    static GLubyte* changeTextureDimension(tgt::ivec3 in_dim, tgt::ivec3 out_dim, GLubyte* indata);
};

/**
 * This class offers an interpolation function for transfer functions. Interpolation: focus on startvalue.
 */
class VRN_CORE_API TransFunc2DPrimitivesStartInterpolationFunction : public TransFunc2DPrimitivesInterpolationFunctionBase {
public:
    TransFunc2DPrimitivesStartInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc2DPrimitivesStartInterpolationFunction"; }
    InterpolationFunction<TransFunc2DPrimitives*>* create() const;
    TransFunc2DPrimitives* interpolate(TransFunc2DPrimitives* startvalue, TransFunc2DPrimitives* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions. Interpolation: focus on endvalue.
 */
class VRN_CORE_API TransFunc2DPrimitivesEndInterpolationFunction : public TransFunc2DPrimitivesInterpolationFunctionBase {
public:
    TransFunc2DPrimitivesEndInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc2DPrimitivesEndInterpolationFunction"; }
    InterpolationFunction<TransFunc2DPrimitives*>* create() const;
    TransFunc2DPrimitives* interpolate(TransFunc2DPrimitives* startvalue, TransFunc2DPrimitives* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions. Interpolation: bisection.
 */
class VRN_CORE_API TransFunc2DPrimitivesStartEndInterpolationFunction : public TransFunc2DPrimitivesInterpolationFunctionBase {
public:
    TransFunc2DPrimitivesStartEndInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc2DPrimitivesStartEndInterpolationFunction"; }
    InterpolationFunction<TransFunc2DPrimitives*>* create() const;
    TransFunc2DPrimitives* interpolate(TransFunc2DPrimitives* startvalue, TransFunc2DPrimitives* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

} // namespace voreen
#endif
