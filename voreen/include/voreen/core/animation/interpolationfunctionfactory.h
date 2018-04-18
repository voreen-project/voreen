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

#ifndef VRN_INTERPOLATIONFUNCTIONFACTORY_H
#define VRN_INTERPOLATIONFUNCTIONFACTORY_H

#include <vector>
#include <string>
#include "voreen/core/voreencoreapi.h"
#include "voreen/core/animation/interpolationfunction.h"
#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

template <class T>
    class InterpolationFunction;

/**
 * Manages all interpolation functions
 */
class VRN_CORE_API InterpolationFunctionFactory : public SerializableFactory {
public:
    InterpolationFunctionFactory();
    ~InterpolationFunctionFactory();

    /**
     * Returns a vector of names of all possible interpolation functions for this type of property.
     */
    template <class T>
    const std::vector<std::string> getListOfFunctionNames() const;

    /**
     * Returns a vector of possible interpolation functions for this type of property.
     */
    template <class T>
    const std::vector<InterpolationFunction<T>*> getListOfFunctions() const;

    /// Returns the default interpolation function for a given property type.
    template <class T>
    const InterpolationFunction<T>* getDefaultFunction() const;
    //TODO!

    /**
     * Returns a new instance of an interpolation function of the given name.
     * If there is no such function, the the default function is returned.
     */
    InterpolationFunctionBase* getFunctionByName(const std::string& name) const;

    /**
     * @see SerializableFactory::getTypeString
     */
    virtual std::string getSerializableTypeString(const std::type_info& type) const;

    /**
     * @see SerializableFactory::createType
     */
    virtual Serializable* createSerializableType(const std::string& className) const;

private:
    /**
     * Registers a new interpolation function to this factory.
     * If the property type is consistent, the given function will be added to the internal list, otherwise it will be deleted.
     */
    void registerFunction(InterpolationFunctionBase* func);

    std::vector<InterpolationFunctionBase*> functions_;

};

template <class T>
const std::vector<std::string> InterpolationFunctionFactory::getListOfFunctionNames() const {
    std::vector<std::string> names;
    typename std::vector<InterpolationFunctionBase*>::const_iterator it;
    for (it = functions_.begin(); it != functions_.end(); ++it) {
        if(dynamic_cast<InterpolationFunction<T>*>(*it))
            names.push_back((*it)->getClassName());
    }

    return names;
}

template <class T>
const std::vector<InterpolationFunction<T>*> InterpolationFunctionFactory::getListOfFunctions() const {
    std::vector<InterpolationFunction<T>*> functions;
    typename std::vector<InterpolationFunctionBase*>::const_iterator it;
    for (it = functions_.begin(); it != functions_.end(); ++it) {
        if(dynamic_cast<InterpolationFunction<T>*>(*it))
            functions.push_back(dynamic_cast<InterpolationFunction<T>*>(*it));
    }
    return functions;
}

template <class T>
const InterpolationFunction<T>* InterpolationFunctionFactory::getDefaultFunction() const {
    return getListOfFunctions<T>().front();
}

} // namespace voreen

#endif
