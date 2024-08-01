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

#include "voreen/core/io/serialization/serializerbase.h"

namespace voreen {

SerializerBase::SerializerBase()
    : errors_()
    , usePointerContentSerialization_(false)
{
}

void SerializerBase::setUsePointerContentSerialization(const bool& usePointerContentSerialization) {
    usePointerContentSerialization_ = usePointerContentSerialization;
}

bool SerializerBase::getUsePointerContentSerialization() const {
    return usePointerContentSerialization_;
}

void SerializerBase::addError(const std::string& message) {
    errors_.push_back(message);
}

void SerializerBase::addError(const std::exception& exception) {
    addError(std::string(exception.what()));
}

void SerializerBase::removeLastError() {
    if (!errors_.empty()) {
        errors_.pop_back();
    }
}

const std::vector<std::string>& SerializerBase::getErrors() const {
    return errors_;
}

SerializerBase::TemporaryUsePointerContentSerializationChanger::TemporaryUsePointerContentSerializationChanger(
    SerializerBase& serializer, const bool &usePointerContentSerialization)
    : serializer_(serializer)
    , storedUsePointerContentSerialization_(serializer_.getUsePointerContentSerialization())
{
    serializer_.setUsePointerContentSerialization(usePointerContentSerialization);
}

SerializerBase::TemporaryUsePointerContentSerializationChanger::~TemporaryUsePointerContentSerializationChanger() {
    serializer_.setUsePointerContentSerialization(storedUsePointerContentSerialization_);
}


} // namespace voreen
