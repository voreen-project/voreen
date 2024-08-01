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

#include "streamline.h"

#include "voreen/core/utils/stringutils.h"

#include <sstream>

namespace voreen {

Streamline::Streamline()
    : magnitudeStatistics_(false)
    , curvatureStatistics_(false)
    , physicalLength_(0.0f)
{
}

Streamline::~Streamline() {
}

//----------------
//  Construction
//----------------
void Streamline::addElementAtEnd(const StreamlineElement& element) {
    if(!streamlineElements_.empty()) {
        // Update physical length.
        physicalLength_ += tgt::distance(streamlineElements_.back().position_, element.position_);

        // Update curvature statistics.
        float currMagnitude = tgt::length(element.velocity_);
        float prevMagnitude = tgt::length(streamlineElements_.back().velocity_);
        if(currMagnitude > 0.0f && prevMagnitude > 0.0f) {
            float angle = std::acos(std::abs(tgt::dot(streamlineElements_.back().velocity_, element.velocity_)) /
                                    (prevMagnitude * currMagnitude));
            curvatureStatistics_.addSample(angle);
        }
    }

    streamlineElements_.push_back(element);
    float length = tgt::length(element.velocity_);
    magnitudeStatistics_.addSample(length);
}

void Streamline::addElementAtFront(const StreamlineElement& element) {
    if(!streamlineElements_.empty()) {
        // Update physical length.
        physicalLength_ += tgt::distance(streamlineElements_.front().position_, element.position_);

        // Update curvature statistics.
        float currMagnitude = tgt::length(element.velocity_);
        float prevMagnitude = tgt::length(streamlineElements_.front().velocity_);
        if(currMagnitude > 0.0f && prevMagnitude > 0.0f) {
            float angle = std::acos(std::abs(tgt::dot(streamlineElements_.front().velocity_, element.velocity_)) /
                                    (prevMagnitude * currMagnitude));
            curvatureStatistics_.addSample(angle);
        }
    }

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
//  Utility
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
            element.position_  = streamlineElements_[pos].position_  * (1.0f - t) + streamlineElements_[pos + 1].position_  * t;
            element.velocity_  = streamlineElements_[pos].velocity_  * (1.0f - t) + streamlineElements_[pos + 1].velocity_  * t;
            element.radius_    = streamlineElements_[pos].radius_    * (1.0f - t) + streamlineElements_[pos + 1].radius_    * t;
            element.time_      = streamlineElements_[pos].time_      * (1.0f - t) + streamlineElements_[pos + 1].time_      * t;

            resampled.addElementAtEnd(element);
        }
    }

    // Append the last element separately.
    resampled.addElementAtEnd(getLastElement());

    return resampled;

}

const Statistics& Streamline::getCurvatureStatistics() const {
    return curvatureStatistics_;
}

const Statistics& Streamline::getMagnitudeStatistics() const {
    return magnitudeStatistics_;
}

float Streamline::getMinMagnitude() const {
    return magnitudeStatistics_.getMin();
}

float Streamline::getMaxMagnitude() const {
    return magnitudeStatistics_.getMax();
}

float Streamline::getPhysicalLength() const {
    return physicalLength_;
}

tgt::vec2 Streamline::getTemporalRange() const {
    return tgt::vec2(getFirstElement().time_, getLastElement().time_);
}

//----------------
//  Storage
//----------------
std::string Streamline::toCSVString(const tgt::mat4& transformationMatrix, const tgt::mat4& velocityTransfomationMatrix) const {
    std::stringstream output;
    output << getNumElements();
    for(const StreamlineElement& element: streamlineElements_) {
        tgt::vec4 transformedPosition = transformationMatrix * tgt::vec4(element.position_,1.f);
        tgt::vec4 transformedVelocity = velocityTransfomationMatrix * tgt::vec4(element.velocity_,1.f);
        output << ", " << transformedPosition.x <<
                  ", " << transformedPosition.y <<
                  ", " << transformedPosition.z <<
                  ", " << transformedVelocity.x <<
                  ", " << transformedVelocity.y <<
                  ", " << transformedVelocity.z <<
                  ", " << element.radius_ <<
                  ", " << element.time_;
    }
    return output.str();
}

void Streamline::serialize(Serializer& s) const {
    //serialize elements as blob
    std::vector<StreamlineElement> vec(streamlineElements_.begin(), streamlineElements_.end());
    s.serializeBinaryBlob("StreamlineElements",vec);
}

void Streamline::deserialize(Deserializer& s) {
    streamlineElements_.clear();

    //deserialize streamlines from binary blob
    std::vector<Streamline::StreamlineElement> vec;
    s.deserializeBinaryBlob("StreamlineElements",vec);
    // We add each element manually since we need to update the statistics.
    for(const StreamlineElement& element : vec) {
        addElementAtEnd(element);
    }
}

}   // namespace
