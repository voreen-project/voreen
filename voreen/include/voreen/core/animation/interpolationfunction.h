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

#ifndef VRN_INTERPOLATIONFUNCTION_H
#define VRN_INTERPOLATIONFUNCTION_H

#include <string>
#include <vector>
#include "voreen/core/animation/interpolationfunctionbase.h"
#include "tgt/camera.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/properties/shaderproperty.h"

namespace voreen {

class Property;
template <class T> class PropertyKeyValue;

/**
 * This class represents an interpolation function between two values.
 */
template <class T>
class InterpolationFunction : public InterpolationFunctionBase {
public:
    InterpolationFunction();
    virtual ~InterpolationFunction();

    virtual InterpolationFunction<T>* create() const = 0;

    /**
     * Interpolates between the given values: start value und endvalue in the time interval [0,1].
     */
    virtual T interpolate(T startvalue, T endvalue, float time) const = 0;

    virtual std::string getGuiName() const = 0;
    virtual std::string getCategory() const = 0;

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);
protected:
    static const std::string loggerCat_;
};

template <class T>
class MultiPointInterpolationFunction : public InterpolationFunction<T> {
public:
    virtual InterpolationFunction<T>* create() const = 0;

    /// Interpolation using several controlpoints.
    virtual T interpolate(std::vector<PropertyKeyValue<T>*> controlpoints, float time) const = 0;
    /// Dummy implementation, use method above.
    virtual T interpolate(T startvalue, T endvalue, float time) const {
        return startvalue;
    }
};

template <class T>
InterpolationFunction<T>::InterpolationFunction() {
}

template <class T>
InterpolationFunction<T>::~InterpolationFunction() {}

template <class T>
void InterpolationFunction<T>::serialize(Serializer& s) const {
}

template <class T>
void InterpolationFunction<T>::deserialize(Deserializer& s) {
}

template <class T>
const std::string InterpolationFunction<T>::loggerCat_("voreen.core.animation.InterpolationFunction");

} // namespace voreen

#endif
