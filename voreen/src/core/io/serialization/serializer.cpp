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

#include "voreen/core/io/serialization/serializer.h"

namespace voreen {

Serializer::Serializer(JsonSerializer& serializer)
    : concreteSerializer(serializer)
{
}
Serializer::Serializer(XmlSerializer& serializer)
    : concreteSerializer(serializer)
{
}

Serializer::~Serializer() {
}

std::string Serializer::getDocumentPath() const {
    DISPATCH_TO_SERIALIZER_RET(getDocumentPath());
}

bool Serializer::getUsePointerContentSerialization() const {
    DISPATCH_TO_SERIALIZER_RET(getUsePointerContentSerialization());
}

void Serializer::setUsePointerContentSerialization(const bool& usePointerContentSerialization) {
    DISPATCH_TO_SERIALIZER(setUsePointerContentSerialization(usePointerContentSerialization));
}

void Serializer::addError(const std::string& message) {
    DISPATCH_TO_SERIALIZER(addError(message));
}

void Serializer::addError(const std::exception& exception) {
    addError(std::string(exception.what()));
}

void Serializer::removeLastError() {
    DISPATCH_TO_SERIALIZER(removeLastError());
}

const std::vector<std::string>& Serializer::getErrors() const {
    DISPATCH_TO_SERIALIZER_RET(getErrors());
}

//Encoding to base64.  Adapted from en.wikibooks.org.
static std::string base64Encode(const std::vector<unsigned char>& inputBuffer) {
    const char encodeLookup[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    const char padCharacter = '=';

    std::string encodedString;
    encodedString.reserve(((inputBuffer.size()/3) + (inputBuffer.size() % 3 > 0)) * 4);
    long temp;
    std::vector<unsigned char>::const_iterator cursor = inputBuffer.begin();
    for(size_t idx = 0; idx < inputBuffer.size()/3; idx++) {
            temp  = (*cursor++) << 16; //Convert to big endian
            temp += (*cursor++) << 8;
            temp += (*cursor++);
            encodedString.append(1,encodeLookup[(temp & 0x00FC0000) >> 18]);
            encodedString.append(1,encodeLookup[(temp & 0x0003F000) >> 12]);
            encodedString.append(1,encodeLookup[(temp & 0x00000FC0) >> 6 ]);
            encodedString.append(1,encodeLookup[(temp & 0x0000003F)      ]);
            //if(encodedString.size() % 72 < 4)
                //encodedString.append(1, '\n');
    }

    switch(inputBuffer.size() % 3) {
        case 1:
                temp  = (*cursor++) << 16; //Convert to big endian
                encodedString.append(1,encodeLookup[(temp & 0x00FC0000) >> 18]);
                encodedString.append(1,encodeLookup[(temp & 0x0003F000) >> 12]);
                encodedString.append(2,padCharacter);
                //if(encodedString.size() % 72 < 3)
                    //encodedString.append(1, '\n');
                break;
        case 2:
                temp  = (*cursor++) << 16; //Convert to big endian
                temp += (*cursor++) << 8;
                encodedString.append(1,encodeLookup[(temp & 0x00FC0000) >> 18]);
                encodedString.append(1,encodeLookup[(temp & 0x0003F000) >> 12]);
                encodedString.append(1,encodeLookup[(temp & 0x00000FC0) >> 6 ]);
                encodedString.append(1,padCharacter);
                //if(encodedString.size() % 72 < 3)
                    //encodedString.append(1, '\n');
                break;
    }
    return encodedString;
}

void Serializer::serializeBinaryBlob(const std::string& key, const unsigned char* data, size_t length) {
    serialize(key, base64Encode(std::vector<unsigned char>(data, data + length)));
}

void Serializer::serializeBinaryBlob(const std::string& key, const std::vector<unsigned char>& data) {
    serialize(key, base64Encode(data));
}


} //namespace voreen
