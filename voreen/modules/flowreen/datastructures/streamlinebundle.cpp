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

#include "streamlinebundle.h"
#include "streamline.h"

#include <sstream>

namespace voreen {

StreamlineBundle::StreamlineBundle()
{
}

StreamlineBundle::StreamlineBundle(Streamline&& prototype)
{
    centroid_ = prototype;
    streamlines_.emplace_back(prototype);
}

void StreamlineBundle::addStreamline(Streamline&& streamline) {

    if(streamline.getNumElements() != centroid_.getNumElements()) {
        tgtAssert(false, "Number of elements does not match");
        return;
    }

    // Update the centroid (using rolling mean).
    Streamline updatedCentroid;
    for(size_t i=0; i<centroid_.getNumElements(); i++) {
        const tgt::vec3& centroidPosition = centroid_.getElementAt(i).position_;
        const tgt::vec3& centroidVelocity = centroid_.getElementAt(i).velocity_;

        tgt::dvec3 deltaPosition = streamline.getElementAt(i).position_ - centroidPosition;
        tgt::dvec3 deltaVelocity = streamline.getElementAt(i).velocity_ - centroidVelocity;

        tgt::vec3 updatedPosition = tgt::dvec3(centroidPosition) + deltaPosition / static_cast<double>(streamlines_.size() + 1);
        tgt::vec3 updatedVelocity = tgt::dvec3(centroidVelocity) + deltaVelocity / static_cast<double>(streamlines_.size() + 1);

        float updatedRadius = std::max(centroid_.getElementAt(i).radius_, 0.5f * tgt::distance(updatedPosition, centroidPosition));

        updatedCentroid.addElementAtEnd(Streamline::StreamlineElement(updatedPosition, updatedVelocity, updatedRadius));
    }
    centroid_ = updatedCentroid;

    // Add the streamline to the bundle.
    streamlines_.emplace_back(streamline);
}

const Streamline& StreamlineBundle::getCentroid() const {
    return centroid_;
}

const std::vector<Streamline>& StreamlineBundle::getStreamlines() const {
    return streamlines_;
}

}
