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

#ifndef VRN_TRANSFUNC1DKEYSINTERPOLATIONFUNCTIONS_H
#define VRN_TRANSFUNC1DKEYSINTERPOLATIONFUNCTIONS_H

#include "voreen/core/animation/interpolationfunction.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API InterpolationFunction<TransFunc1DKeys*>;
#endif

/// Base class for transfer function interpolation, provides commonly used methods.
class VRN_CORE_API TransFunc1DKeysInterpolationFunctionBase : public InterpolationFunction<TransFunc1DKeys*> {
protected:
    static GLubyte* convertTextureToRGBA(tgt::ivec3 dim, GLubyte* textur, GLuint inputformat);
    static GLubyte* changeTextureDimension(tgt::ivec3 in_dim, tgt::ivec3 out_dim, GLubyte* indata);
};

/// Default interpolation
class VRN_CORE_API TransFunc1DKeysInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions. Interpolation: focus on startvalue.
 */
class VRN_CORE_API TransFunc1DKeysStartInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysStartInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysStartInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions. Interpolation: focus on endvalue.
 */
class VRN_CORE_API TransFunc1DKeysEndInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysEndInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysEndInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions. Interpolation: bisection.
 */
class VRN_CORE_API TransFunc1DKeysStartEndInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysStartEndInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysStartEndInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys
 * the functions will be interpolated linearly (keywise).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quadratically (keywise, easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuadInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuadInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuadInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quadratically (keywise, easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuadOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuadOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuadOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quadratically (keywise, easing in, then easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuadInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuadInOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuadInOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quadratically (keywise, easing out, then easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuadOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuadOutInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuadOutInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated cubicularly (keywise, easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseCubicInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseCubicInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseCubicInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated cubicularly (keywise, easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseCubicOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseCubicOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseCubicOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated cubicularly (keywise, easing in, then easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseCubicInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseCubicInOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseCubicInOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated cubicularly (keywise, easing out, then easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseCubicOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseCubicOutInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseCubicOutInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quartetically (keywise, easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuartInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuartInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuartInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quartetically (keywise, easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuartOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuartOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuartOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quartetically (keywise, easing in, then easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuartInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuartInOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuartInOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quartetically (keywise, easing out, then easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuartOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuartOutInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuartOutInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quintically (keywise, easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuintInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuintInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuintInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quintically (keywise, easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuintOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuintOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuintOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quintically (keywise, easing in, then easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuintInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuintInOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuintInOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated quintically (keywise, easing out, then easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseQuintOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseQuintOutInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseQuintOutInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated sineousidally (keywise, easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseSineInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseSineInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseSineInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated sineousidally (keywise, easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseSineOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseSineOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseSineOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated sineousidally (keywise, easing in, then easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseSineInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseSineInOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseSineInOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated sineousidally (keywise, easing out, then easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseSineOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseSineOutInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseSineOutInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated exponentially (keywise, easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseExponentInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseExponentInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseExponentInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated exponentially (keywise, easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseExponentOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseExponentOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseExponentOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated exponentially (keywise, easing in, then easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseExponentInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseExponentInOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseExponentInOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated exponentially (keywise, easing out, then easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseExponentOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseExponentOutInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseExponentOutInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated circularly (keywise, easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseCircInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseCircInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseCircInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated circularly (keywise, easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseCircOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseCircOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseCircOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated circularly (keywise, easing in, then easing out).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseCircInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseCircInOutInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseCircInOutInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

/**
 * This class offers an interpolation function for transfer functions.
 * If the startfunction and the endfunction are both 1D-functions with the same number of keys,
 * the functions will be interpolated circularly (keywise, easing out, then easing in).
 * If not this functions falls back to a default function like TransFunc1DKeysStartInterpolationFunction.
 */
class VRN_CORE_API TransFunc1DKeysKeyWiseCircOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
public:
    TransFunc1DKeysKeyWiseCircOutInInterpolationFunction();
    virtual std::string getClassName() const { return "TransFunc1DKeysKeyWiseCircOutInInterpolationFunction"; }
    InterpolationFunction<TransFunc1DKeys*>* create() const;
    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;

    std::string getGuiName() const;
    std::string getCategory() const;
};

///**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated linearly.
// */
//class VRN_CORE_API TransFuncTextureLinearInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureLinearInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureLinearInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quadratically (easing-in).
// */
//class VRN_CORE_API TransFuncTextureQuadInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuadInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuadInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quadratically (easing-out).
// */
//class VRN_CORE_API TransFuncTextureQuadOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuadOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuadOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quadratically (easing-in, then easing-out).
// */
//class VRN_CORE_API TransFuncTextureQuadInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuadInOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuadInOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quadratically (easing-out, then easing-in).
// */
//class VRN_CORE_API TransFuncTextureQuadOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuadOutInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuadOutInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated cubicularly (easing-in).
// */
//class VRN_CORE_API TransFuncTextureCubicInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureCubicInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureCubicInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated cubicularly (easing-out).
// */
//class VRN_CORE_API TransFuncTextureCubicOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureCubicOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureCubicOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated cubicularly (easing-in, then easing-out).
// */
//class VRN_CORE_API TransFuncTextureCubicInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureCubicInOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureCubicInOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated cubicularly (easing-out, then easing-in).
// */
//class VRN_CORE_API TransFuncTextureCubicOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureCubicOutInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureCubicOutInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quartetically (easing-in).
// */
//class VRN_CORE_API TransFuncTextureQuartInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuartInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuartInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quartetically (easing-out).
// */
//class VRN_CORE_API TransFuncTextureQuartOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuartOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuartOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quartetically (easing-in, then easing-out).
// */
//class VRN_CORE_API TransFuncTextureQuartInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuartInOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuartInOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quartetically (easing-out, then easing-in).
// */
//class VRN_CORE_API TransFuncTextureQuartOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuartOutInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuartOutInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quintically (easing-in).
// */
//class VRN_CORE_API TransFuncTextureQuintInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuintInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuintInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quintically (easing-out).
// */
//class VRN_CORE_API TransFuncTextureQuintOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuintOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuintOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quintically (easing-in, then easing-out).
// */
//class VRN_CORE_API TransFuncTextureQuintInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuintInOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuintInOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated quintically (easing-out, then easing-in).
// */
//class VRN_CORE_API TransFuncTextureQuintOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureQuintOutInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureQuintOutInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated sineousidally (easing-in).
// */
//class VRN_CORE_API TransFuncTextureSineInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureSineInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureSineInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated sineousidally (easing-out).
// */
//class VRN_CORE_API TransFuncTextureSineOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureSineOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureSineOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated sineousidally (easing-in, then easing-out).
// */
//class VRN_CORE_API TransFuncTextureSineInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureSineInOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureSineInOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated sineousidally (easing-out, then easing-in).
// */
//class VRN_CORE_API TransFuncTextureSineOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureSineOutInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureSineOutInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated exponentially (easing-in).
// */
//class VRN_CORE_API TransFuncTextureExponentInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureExponentInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureExponentInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated exponentially (easing-out).
// */
//class VRN_CORE_API TransFuncTextureExponentOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureExponentOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureExponentOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated exponentially (easing-in, then easing-out).
// */
//class VRN_CORE_API TransFuncTextureExponentInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureExponentInOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureExponentInOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated exponentially (easing-out, then easing-in).
// */
//class VRN_CORE_API TransFuncTextureExponentOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureExponentOutInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureExponentOutInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated circularly (easing-in).
// */
//class VRN_CORE_API TransFuncTextureCircInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureCircInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureCircInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated circularly (easing-out).
// */
//class VRN_CORE_API TransFuncTextureCircOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureCircOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureCircOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated circularly (easing-in, then easing-out).
// */
//class VRN_CORE_API TransFuncTextureCircInOutInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureCircInOutInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureCircInOutInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};
//
//**
// * This class offers an interpolation function for transfer functions.
// * The two textures of the transferfunctions are interpolated circularly (easing-out, then easing-in).
// */
//class VRN_CORE_API TransFuncTextureCircOutInInterpolationFunction : public TransFunc1DKeysInterpolationFunctionBase {
//public:
//    TransFuncTextureCircOutInInterpolationFunction();
//    virtual std::string getClassName() const { return "TransFuncTextureCircOutInInterpolationFunction"; }
//    InterpolationFunction<TransFunc1DKeys*>* create() const;
//    TransFunc1DKeys* interpolate(TransFunc1DKeys* startvalue, TransFunc1DKeys* endvalue, float time) const;
//
//    std::string getGuiName() const;
//    std::string getCategory() const;
//};

} // namespace voreen
#endif
