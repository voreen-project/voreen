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

#ifndef VRN_ANNOTATION_H
#define VRN_ANNOTATION_H

#include "voreen/core/utils/GLSLparser/token.h"

#include <sstream>
#include <string>
#include <vector>

namespace voreen {

namespace glslparser {

class GLSLAnnotation {
public:
    GLSLAnnotation(const std::string& name, const std::vector<Token*>& value)
        : name_(name),
        value_(value.size(), 0)
    {
        for (size_t i = 0; i < value.size(); ++i)
            value_[i] = value[i]->getCopy();
    }

    virtual ~GLSLAnnotation() {
        for (size_t i = 0; i < value_.size(); ++i)
            delete value_[i];
    }

    const std::string& getName() const { return name_; }
    size_t getNumValueElements() { return value_.size(); }

    template<typename T>
    std::vector<T> getValueAs() const {
        std::vector<T> ret;
        for (size_t i = 0; i < value_.size(); ++i) {
            GenericToken<std::string>* const gen =
                dynamic_cast<GenericToken<std::string>* const>(value_[i]);

            if (gen == 0)
                continue;

            std::istringstream iss(gen->getValue());
            T aux(0);
            iss >> aux;
            ret.push_back(aux);
        }

        return ret;
    }

    std::string getValueString() const {
        std::ostringstream oss;
        for (size_t i = 0; i < value_.size(); ++i) {
            GenericToken<std::string>* const gen =
                dynamic_cast<GenericToken<std::string>* const>(value_[i]);

            if (gen == 0)
                continue;

            if (i > 0)
                oss << ", ";
            oss << gen->getValue();
        }
        return oss.str();
    }

protected:
    std::string name_;
    std::vector<Token*> value_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
