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

#include "cmmergertree.h"

namespace voreen {


//CMMergerTree
const int CMMergerTree::NO_HALO_ID = -1;

CMMergerTree::CMMergerTree()
    : containsGalacticusData_(false)
{
}

CMMergerTree::~CMMergerTree()
{
}

bool CMMergerTree::containsGalacticusData() const {
    return containsGalacticusData_;
}
const CMTimeStepHaloList* CMMergerTree::haloDataAt(float a) const {
    if(stepListsBegin() == stepListsEnd()) {
        return nullptr;
    }
    const CMTimeStepHaloList* lastData = stepListsBegin();
    //We know that HaloList_ is sorted, so we can just go from left to right:
    for(const CMTimeStepHaloList* currentData = lastData+1; currentData != stepListsEnd(); ++currentData) {
        if(lastData->a + currentData->a > 2*a) {
            if(lastData == nullptr) {
                return currentData;
            } else {
                return lastData;
            }
        } else {
            lastData = currentData;
        }
    }
    return lastData;
}

const CMHalo* CMMergerTree::haloByID(int id) const {
    const CMHalo* address = halosBegin()+id;
    if(id < 0 || address>=halosEnd()) {
        return nullptr;
    }
    return address;
}


//CMTimeStepHaloList
const CMTimeStepHaloList* CMMergerTree::stepListsBegin() const {
    return reinterpret_cast<const CMTimeStepHaloList*>(reinterpret_cast<const char*>(this)+hlistBeginOffset);
}
const CMTimeStepHaloList* CMMergerTree::stepListsEnd() const {
    return reinterpret_cast<const CMTimeStepHaloList*>(reinterpret_cast<const char*>(this)+hlistEndOffset);
}
const CMHalo* CMMergerTree::halosBegin() const {
    return reinterpret_cast<const CMHalo*>(reinterpret_cast<const char*>(this)+halosBeginOffset);
}
const CMHalo* CMMergerTree::halosEnd() const {
    return reinterpret_cast<const CMHalo*>(reinterpret_cast<const char*>(this)+halosEndOffset);
}
const CMHalo* CMTimeStepHaloList::halosBegin() const {
    return reinterpret_cast<const CMHalo*>(reinterpret_cast<const char*>(this)+beginOffset);
}
const CMHalo* CMTimeStepHaloList::halosEnd() const {
    return reinterpret_cast<const CMHalo*>(reinterpret_cast<const char*>(this)+endOffset);
}
size_t CMTimeStepHaloList::size() const {
    return std::distance(halosBegin(), halosEnd());
}


//CMHalo
const CMHalo* CMHalo::unsafehaloByID(int id) const {
    if(id==CMMergerTree::NO_HALO_ID) {
        return nullptr;
    }
    return this-this->ID+id;
}
const CMHalo* CMHalo::descendant() const {
    return unsafehaloByID(descendantID);
}
const CMHalo* CMHalo::host() const {
    return unsafehaloByID(hostID);
}
const CMHalo* CMHalo::rootHost() const {
    return unsafehaloByID(rootHostID);
};
const CMHalo* CMHalo::parent() const {
    return unsafehaloByID(parentID);
}
const CMHalo* CMHalo::spouse() const {
    return unsafehaloByID(spouseID);
}
const CMHalo* CMHalo::satellite() const {
    return unsafehaloByID(satelliteID);
}
const CMHalo* CMHalo::siblingSatellite() const {
    return unsafehaloByID(siblingSatelliteID);
}
}
