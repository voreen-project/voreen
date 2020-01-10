/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "streamline.h"

#include "voreen/core/utils/stringutils.h"

#include "tgt/assert.h"

#include <sstream>

namespace voreen {

Streamline::Streamline()
    : magnitudeStatistics_(false)
{
}

Streamline::~Streamline() {
}

//----------------
//  Construction
//----------------
void Streamline::addElementAtEnd(const StreamlineElement& element) {
    streamlineElements_.push_back(element);
    float length = tgt::length(element.velocity_);
    magnitudeStatistics_.addSample(length);
}

void Streamline::addElementAtFront(const StreamlineElement& element) {
    streamlineElements_.push_front(element);
    float length = tgt::length(element.velocity_);
    magnitudeStatistics_.addSample(length);
}

//----------------
//  Access
//----------------
size_t Streamline::getNumElements() const {
    return streamlineElements_.size();
}

float Streamline::getMinMagnitude() const {
    return magnitudeStatistics_.getMin();
}

float Streamline::getMaxMagnitude() const {
    return magnitudeStatistics_.getMax();
}

const Streamline::StreamlineElement& Streamline::getElementAt(size_t pos) const {
    tgtAssert(pos < getNumElements(), "Requested element does not exist!");
    return streamlineElements_[pos];
}

const Streamline::StreamlineElement& Streamline::getFirstElement() const {
    return streamlineElements_.front();
}

const Streamline::StreamlineElement& Streamline::getLastElement() const {
    return streamlineElements_.back();
}

//----------------
//  Utilitiy
//----------------
Streamline Streamline::resample(size_t samples) const {

    // This will hold the resampled streamline.
    Streamline resampled;

    // Append the first element separately.
    resampled.addElementAtEnd(getFirstElement());

    // In case we have more than two samples, we achieve the result by linear interpolation.
    if(samples > 2) {

        std::vector<float> distances(getNumElements());

        float totalLength = 0.0f;
        distances[0] = 0.0f;
        for(size_t i = 1; i < getNumElements(); i++) {
            totalLength += tgt::distance(streamlineElements_[i-1].position_, streamlineElements_[i].position_);
            distances[i] = totalLength;
        }

        const float segmentLength = totalLength / (samples - 1);
        for (size_t i = 1, pos = 0; i < samples - 1; i++) {

            while (distances[pos+1] < segmentLength * i) pos++;

            const float t = tgt::clamp((segmentLength * i - distances[pos]) / (distances[pos + 1] - distances[pos]), 0.0f, 1.0f);

            Streamline::StreamlineElement element;
            element.position_ = streamlineElements_[pos].position_ * (1.0f - t) + streamlineElements_[pos + 1].position_ * t;
            element.velocity_ = streamlineElements_[pos].velocity_ * (1.0f - t) + streamlineElements_[pos + 1].velocity_ * t;

            resampled.addElementAtEnd(element);
        }
    }

    // Append the last element separately.
    resampled.addElementAtEnd(getLastElement());

    return resampled;

}

//----------------
//  Storage
//----------------
std::string Streamline::toCSVString(const tgt::mat4& transformationMatrix, const tgt::mat4& velocityTransfomationMatrix) const {
    std::stringstream output;
    output << getNumElements() << ", " << getMinMagnitude() << ", " << getMaxMagnitude();
    for(size_t i = 0; i < streamlineElements_.size(); i++) {
        tgt::vec4 transformedPosition = transformationMatrix * tgt::vec4(streamlineElements_[i].position_,1.f);
        tgt::vec4 transformedVelocity = velocityTransfomationMatrix * tgt::vec4(streamlineElements_[i].velocity_,1.f);
        output << ", " << transformedPosition.x <<
                  ", " << transformedPosition.y <<
                  ", " << transformedPosition.z <<
                  ", " << transformedVelocity.x <<
                  ", " << transformedVelocity.y <<
                  ", " << transformedVelocity.z;
    }
    return output.str();
}

void Streamline::serialize(Serializer& s) const {
    s.serialize("minMagnitude_",magnitudeStatistics_.getMin());
    s.serialize("maxMagnitude_",magnitudeStatistics_.getMax());
    //serialize elements as blob
    std::vector<StreamlineElement> vec(streamlineElements_.size());
    std::copy(streamlineElements_.begin(),streamlineElements_.end(),vec.begin());
    s.serializeBinaryBlob("StreamlineElements",vec);
}

void Streamline::deserialize(Deserializer& s) {
    //s.deserialize("minMagnitude_",minMagnitude_);
    //s.deserialize("maxMagnitude_",maxMagnitude_);
    //deserialize streamlines from binary blob
    std::vector<Streamline::StreamlineElement> vec;
    s.deserializeBinaryBlob("StreamlineElements",vec);
    // We add each element manually since we need to update the statistics.
    streamlineElements_.clear();
    for(const StreamlineElement& element : vec) {
        addElementAtEnd(element);
    }
}

}   // namespace
