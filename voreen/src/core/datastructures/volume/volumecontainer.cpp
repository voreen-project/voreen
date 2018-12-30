/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/core/datastructures/volume/volumecontainer.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/io/progressbar.h"

namespace voreen {

const std::string VolumeContainer::loggerCat_ = "voreen.VolumeContainer";

VolumeContainer::VolumeContainer() :
    VolumeList()
{}

VolumeContainer::~VolumeContainer() {
    //clear(); // Clear is called in super class destructor.
}

void VolumeContainer::remove(const VolumeBase* handle) {
    if (contains(handle)) {
        VolumeList::remove(handle); // also notifies the observers
        delete handle;
    }
}

void VolumeContainer::remove(const VolumeList* volumeList) {
    tgtAssert(volumeList, "Unexpected null pointer");

    // the passed parameter might be object itself,
    // therefore we have to copy over the volumes to a
    // temporary vector
    std::vector<VolumeBase*> deleteList;
    for (size_t i=0; i<volumeList->size(); ++i) {
        deleteList.push_back(volumeList->at(i));
    }

    for (size_t i=0; i<deleteList.size(); ++i) {
        remove(deleteList[i]); // also notifies the observers
    }
}

void VolumeContainer::clear() {
    while (!empty())
        remove(first());
}

void VolumeContainer::release(const VolumeBase* handle) {
    VolumeList::remove(handle);
}

void VolumeContainer::releaseAll() {
    VolumeList::clear();
}

} // namespace
