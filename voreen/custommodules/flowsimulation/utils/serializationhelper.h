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

#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/io/serialization/deserializer.h"

#include "../ext/openlb/voreen/openlb_parameters.h"

namespace voreen {

template<typename T, typename S>
void serializeVector(Serializer& s, const std::string& key, const std::vector<S>& input) {
    std::vector<T> result(input.size());
    for(size_t i=0; i<input.size(); i++) {
        result[i] = static_cast<T>(input[i]); // Enforced slicing!
    }

    s.template serialize(key, result);
}

template<typename T, typename S>
void deserializeVector(Deserializer& s, const std::string& key, std::vector<S>& output) {
    std::vector<T> result;
    s.template deserialize(key, result);

    output.resize(result.size());
    for(size_t i=0; i<output.size(); i++) {
        output[i] = static_cast<S>(result[i]); // Enforced slicing!
    }
}

class VelocityCurveSerializable : public VelocityCurve, public Serializable {
public:

    VelocityCurveSerializable();
    VelocityCurveSerializable(const VelocityCurve& other);

    void serialize(Serializer& s) const;
    void deserialize(Deserializer& s);
};

class FlowIndicatorSerializable : public FlowIndicator, public Serializable {
public:

    FlowIndicatorSerializable();
    FlowIndicatorSerializable(const FlowIndicator& other);

    void serialize(Serializer& s) const;
    void deserialize(Deserializer& s);
};

class ParametersSerializable : public Parameters, public Serializable {
public:
    ParametersSerializable();
    ParametersSerializable(const Parameters& other);

    void serialize(Serializer& s) const;
    void deserialize(Deserializer& s);
};

}