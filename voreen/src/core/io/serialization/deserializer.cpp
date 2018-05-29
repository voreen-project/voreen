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

#include "voreen/core/io/serialization/deserializer.h"

namespace voreen {

Deserializer::Deserializer(XmlDeserializer& deserializer)
    : concreteDeserializer(deserializer)
{
}
Deserializer::Deserializer(JsonDeserializer& deserializer)
    : concreteDeserializer(deserializer)
{
}

Deserializer::~Deserializer() {
}

std::string Deserializer::getDocumentPath() const {
    DISPATCH_TO_DESERIALIZER_RET(getDocumentPath());
}

bool Deserializer::getUsePointerContentSerialization() const {
    DISPATCH_TO_DESERIALIZER_RET(getUsePointerContentSerialization());
}

void Deserializer::setUsePointerContentSerialization(const bool& usePointerContentSerialization) {
    DISPATCH_TO_DESERIALIZER(setUsePointerContentSerialization(usePointerContentSerialization));
}

void Deserializer::addError(const std::string& message) {
    DISPATCH_TO_DESERIALIZER(addError(message));
}

void Deserializer::addError(const std::exception& exception) {
    addError(std::string(exception.what()));
}

void Deserializer::removeLastError() {
    DISPATCH_TO_DESERIALIZER(removeLastError());
}

const std::vector<std::string>& Deserializer::getErrors() const {
    DISPATCH_TO_DESERIALIZER_RET(getErrors());
}

static std::vector<unsigned char> base64Decode(const std::string& input) {
    //input.erase(std::remove(input.begin(), input.end(), '\n'), input.end());

    const char padCharacter = '=';

    if (input.length() % 4) //Sanity check
        throw SerializationFormatException("Non-Valid base64!");

    size_t padding = 0;
    if (input.length()) {
        if (input[input.length()-1] == padCharacter)
                padding++;
        if (input[input.length()-2] == padCharacter)
                padding++;
    }

    //Setup a vector to hold the result
    std::vector<unsigned char> decodedBytes;
    decodedBytes.reserve(((input.length()/4)*3) - padding);
    long temp = 0; //Holds decoded quanta
    std::string::const_iterator cursor = input.begin();
    while (cursor < input.end()) {
        for (size_t quantumPosition = 0; quantumPosition < 4; quantumPosition++) {
            temp <<= 6;
            if       (*cursor >= 0x41 && *cursor <= 0x5A) // This area will need tweaking if
                temp |= *cursor - 0x41;                       // you are using an alternate alphabet
            else if  (*cursor >= 0x61 && *cursor <= 0x7A)
                temp |= *cursor - 0x47;
            else if  (*cursor >= 0x30 && *cursor <= 0x39)
                temp |= *cursor + 0x04;
            else if  (*cursor == 0x2B)
                temp |= 0x3E; //change to 0x2D for URL alphabet
            else if  (*cursor == 0x2F)
                temp |= 0x3F; //change to 0x5F for URL alphabet
            else if  (*cursor == padCharacter) { //pad
                switch( input.end() - cursor ) {
                    case 1: //One pad character
                        decodedBytes.push_back((temp >> 16) & 0x000000FF);
                        decodedBytes.push_back((temp >> 8 ) & 0x000000FF);
                        return decodedBytes;
                    case 2: //Two pad characters
                        decodedBytes.push_back((temp >> 10) & 0x000000FF);
                        return decodedBytes;
                    default:
                        throw SerializationFormatException("Invalid Padding in Base 64!");
                }
            }  else
                throw SerializationFormatException("Non-Valid Character in Base 64!");
            cursor++;
        }
        decodedBytes.push_back((temp >> 16) & 0x000000FF);
        decodedBytes.push_back((temp >> 8 ) & 0x000000FF);
        decodedBytes.push_back((temp      ) & 0x000000FF);
    }

    return decodedBytes;
}

void Deserializer::deserializeBinaryBlob(const std::string& key, unsigned char* inputBuffer, size_t reservedMemory) {
    std::string tmp;
    deserialize(key, tmp);
    std::vector<unsigned char> dataVec = base64Decode(tmp);

    if(dataVec.size() != reservedMemory)
        throw SerializationFormatException("Json node with key '" + key + "' contains a binary blob for which unsufficient memory has been reserved.");
    else {
        //std::copy(inputBuffer, inputBuffer + dataVec.size(), dataVec.begin()); //doesn't work
        memcpy(inputBuffer, &dataVec[0], reservedMemory);
    }
}

void Deserializer::deserializeBinaryBlob(const std::string& key, std::vector<unsigned char>& buffer) {
    std::string tmp;
    deserialize(key, tmp);
    buffer = base64Decode(tmp);
}

} // namespace voreen
