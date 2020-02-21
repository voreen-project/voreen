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

#include "voreen/core/utils/voreenqualitymode.h"

namespace voreen {

VoreenQualityMode::VoreenQualityMode() : quality_(RQ_DEFAULT){}

VoreenQualityMode::~VoreenQualityMode() {}

    //----------------------------
    //        Getter
    //----------------------------
VoreenQualityMode::RenderingQuality VoreenQualityMode::getQuality() const {
    return quality_;
}

bool VoreenQualityMode::isInteractionMode() const {
    return (quality_ == RQ_INTERACTIVE);
}

bool VoreenQualityMode::isDefaultMode() const {
    return (quality_ == RQ_DEFAULT);
}

bool VoreenQualityMode::isHighQualityMode() const {
    return (quality_ == RQ_HIGH);
}

void VoreenQualityMode::setQuality(RenderingQuality quality) {
    if(quality_ != quality) {
        quality_ = quality;
        notifyQualiyModeChanged();
        LDEBUGC("voreen.VoreenQualityMode", "Quality mode changed to: " << quality);
    }
}

    //----------------------------
    //    Requests
    //----------------------------
void VoreenQualityMode::requestQualityMode(RenderingQuality requestedQuality, void* source) {
    switch (requestedQuality) {
    case RQ_DEFAULT:
        //current mode is not longer needed
        activeSources_.erase(source);
        if (activeSources_.empty())
            setQuality(RQ_DEFAULT);
        break;
    case RQ_INTERACTIVE:
    case RQ_HIGH:
        if(quality_ == requestedQuality) {
            //add source
            activeSources_.insert(source);
        } else {
            //switch quality if it was default
            if(quality_ == RQ_DEFAULT) {
                activeSources_.insert(source);
                setQuality(requestedQuality);
            } else {
                //ignore request
            }
        }
        break;
    default:
        tgtAssert(false,"Unknown quality!");
    }
}

void VoreenQualityMode::requestQualityModeForced(RenderingQuality requestedQuality, void* source) {
    switch (requestedQuality) {
    case RQ_DEFAULT:
        //current mode is not longer needed
        activeSources_.erase(source);
        if (activeSources_.empty())
            setQuality(RQ_DEFAULT);
        break;
    case RQ_INTERACTIVE:
    case RQ_HIGH:
        if(quality_ == requestedQuality) {
            //add source
            activeSources_.insert(source);
        } else {
            //switch quality if it was default
            if(quality_ == RQ_DEFAULT) {
                activeSources_.insert(source);
                setQuality(requestedQuality);
            } else {
                //force request
                activeSources_.clear();
                activeSources_.insert(source);
                setQuality(requestedQuality);
            }
        }
        break;
    default:
        tgtAssert(false,"Unknown quality!");
    }
}

void VoreenQualityMode::notifyQualiyModeChanged() {
    const std::vector<QualityModeObserver*> vec = getObservers();
    for(size_t i = 0; i < vec.size(); i++)
        vec[i]->qualityModeChanged();
}

} // namespace
