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

#include "streamlinebundle.h"
#include "streamline.h"

#include "voreen/core/io/serialization/serialization.h"

#include <sstream>

namespace voreen {

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
//      StreamlineBundle
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
StreamlineBundle::StreamlineBundle()
    : radius_(0.0f)
{
}

StreamlineBundle::StreamlineBundle(size_t index, const Streamline& prototype)
    : radius_(0.0f)
{
    nodes_.reserve(prototype.getNumElements());
    for(size_t i = 0; i < prototype.getNumElements(); i++) {
        Node node;
        node.position_ = prototype.getElementAt(i).position_;
        node.velocity_ = prototype.getElementAt(i).velocity_;
        nodes_.push_back(node);
    }

    streamlines_.push_back(index);
}

void StreamlineBundle::addStreamline(size_t index, const Streamline& streamline) {

    // Get current centroid.
    Streamline centroid = getCentroid();

    // Update the approximate radius.
    float minimalDistanceSq = std::numeric_limits<float>::max();
    for(size_t k = 0; k < nodes_.size(); k++) {
        for(size_t j = 0; j < streamline.getNumElements(); j++) {
            minimalDistanceSq = std::min(minimalDistanceSq, tgt::distanceSq(centroid.getElementAt(k).position_, streamline.getElementAt(j).position_));
        }

        nodes_[k].position_ += streamline.getElementAt(k).position_;
        nodes_[k].velocity_ += streamline.getElementAt(k).velocity_;
    }

    // Update radius.
    radius_ += sqrtf(minimalDistanceSq);

    // Add streamline to bundle
    streamlines_.push_back(index);
}

const std::vector<size_t>& StreamlineBundle::getStreamlines() const {
    return streamlines_;
}

Streamline StreamlineBundle::getCentroid() const {

    Streamline streamline;

    for(size_t i = 0; i < nodes_.size(); i++) {
        Streamline::StreamlineElement element;
        element.position_ = nodes_[i].position_ / static_cast<float>(streamlines_.size());
        element.velocity_ = nodes_[i].velocity_ / static_cast<float>(streamlines_.size());
        streamline.addElementAtEnd(element);
    }

    return streamline;
}

float StreamlineBundle::getRadius() const {
    return radius_ / streamlines_.size();
}

//----------------
//  Storage
//----------------
std::string StreamlineBundle::toCSVString(const tgt::mat4& transformationMatrix) const {
    std::stringstream output;

    output << streamlines_.size() << ", " << getRadius() << ", ";
    output << getCentroid().toCSVString(transformationMatrix);

    return output.str();
}

void StreamlineBundle::serialize(Serializer& s) const {
    s.serialize("Radius", radius_);
    s.serializeBinaryBlob("AbsorbedStreamlines", streamlines_);
    s.serializeBinaryBlob("StreamlineBundleNodes", nodes_);
}

void StreamlineBundle::deserialize(Deserializer& s) {
    s.deserialize("Radius", radius_);
    s.deserializeBinaryBlob("AbsorbedStreamlines", streamlines_);
    s.deserializeBinaryBlob("StreamlineBundleNodes", nodes_);
}

}
