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

#include "voreen/core/datastructures/volume/volumelist.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/utils/hashing.h"

using std::vector;

namespace voreen {

const std::string VolumeList::loggerCat_ = "voreen.VolumeList";

VolumeList::VolumeList()
    : Observable<VolumeListObserver>()
{
}

VolumeList::~VolumeList() {
    clear();
}

void VolumeList::add(VolumeBase* volume) {

    tgtAssert(volume, "Null pointer as Volume passed");
    if (!contains(volume)) {
        notifyPendingDataInvalidation();
        volumes_.push_back(volume);
        volume->Observable<VolumeObserver>::addObserver(static_cast<const VolumeObserver*>(this));
        volume->DataInvalidationObservable::addObserver(static_cast<const DataInvalidationObserver*>(this));
        notifyVolumeAdded(volume);
    }
}

void VolumeList::add(const VolumeList* volumeList) {

    tgtAssert(volumeList, "Unexpected null pointer");
    for (size_t i=0; i<volumeList->size(); ++i) {
        add(volumeList->at(i));
    }
}

void VolumeList::remove(const VolumeBase* volume){
    std::vector<VolumeBase*>::iterator handleIter = find(volume);
    if (handleIter != volumes_.end()) {
        notifyPendingDataInvalidation();
        volumes_.erase(handleIter);
        volume->DataInvalidationObservable::removeObserver(static_cast<const DataInvalidationObserver*>(this));
        volume->Observable<VolumeObserver>::removeObserver(static_cast<const VolumeObserver*>(this));
        notifyVolumeRemoved(volume);
    }
}

void VolumeList::remove(const VolumeList* volumeList) {
   tgtAssert(volumeList, "Unexpected null pointer");
    for (size_t i=0; i<volumeList->size(); ++i) {
        remove(volumeList->at(i));
    }
}

std::vector<VolumeBase*>::iterator VolumeList::find(const VolumeBase* volume) {
    return std::find(volumes_.begin(), volumes_.end(), volume);
}

std::vector<VolumeBase*>::const_iterator VolumeList::find(const VolumeBase* volume) const {
    return std::find(volumes_.begin(), volumes_.end(), volume);
}

bool VolumeList::contains(const VolumeBase* volume) const {
    return (find(volume) != volumes_.end());
}

VolumeBase* VolumeList::at(size_t i) const {
    tgtAssert(i < volumes_.size(), "Invalid index");
    return volumes_.at(i);
}

VolumeBase* VolumeList::first() const {
    return (!empty() ? volumes_.front() : 0);
}

void VolumeList::clear() {
    while (!empty())
        remove(first());
}

VolumeList* VolumeList::selectModality(const Modality& modality) const {

    VolumeList* collection = new VolumeList();
    for (size_t i=0; i<volumes_.size(); ++i) {
        if (volumes_[i]->getModality() == modality)
            collection->add(volumes_[i]);
    }
    return collection;
}

voreen::VolumeList* VolumeList::selectOrigin(const VolumeURL& origin) const {
    VolumeList* collection = new VolumeList();
    for (size_t i=0; i<volumes_.size(); ++i) {
        Volume* vh = dynamic_cast<Volume*>(volumes_[i]);
        if (vh && vh->getOrigin() == origin)
            collection->add(volumes_[i]);
    }
    return collection;
}

VolumeList* VolumeList::subList(size_t start, size_t end) const {
    VolumeList* subCollection = new VolumeList();
    tgtAssert(start <= end, "invalid indices");
    tgtAssert(start < volumes_.size(), "invalid start index");
    tgtAssert(end < volumes_.size(), "invalid end index");
    for (size_t index = start; index <= end; index++)
        subCollection->add(volumes_.at(index));
    return subCollection;
}

VolumeList* VolumeList::subList(const std::vector<size_t>& indices) const {
    VolumeList* subCollection = new VolumeList();
    for (size_t i=0; i<indices.size(); i++) {
        tgtAssert(indices.at(i) < volumes_.size(), "invalid index");
        subCollection->add(volumes_.at(indices.at(i)));
    }
    return subCollection;
}


size_t VolumeList::size() const {
    return volumes_.size();
}

bool VolumeList::empty() const {
    return (size() == 0);
}

std::string VolumeList::getHash() const {
    if(empty())
        return VoreenHash::getHash("emptylist");
    std::stringstream stream;
    for( std::vector<VolumeBase*>::const_iterator vol = volumes_.begin(); vol != volumes_.end(); vol++)
        stream << (*vol)->getHash();
    return VoreenHash::getHash(stream.str());
}

void VolumeList::notifyVolumeAdded(const VolumeBase* handle) {
    const vector<VolumeListObserver*> observers = Observable<VolumeListObserver>::getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->volumeAdded(this, handle);

}

void VolumeList::notifyVolumeRemoved(const VolumeBase* handle) {
    const vector<VolumeListObserver*> observers = Observable<VolumeListObserver>::getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->volumeRemoved(this, handle);
}

void VolumeList::notifyVolumeChanged(const VolumeBase* handle) {
    const vector<VolumeListObserver*> observers = Observable<VolumeListObserver>::getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->volumeChanged(this, handle);
}

// implementation of VolumeObserver interface
void VolumeList::volumeChange(const VolumeBase* handle) {
    if (contains(handle))
        notifyVolumeChanged(handle);

}

// implementation of VolumeObserver interface
void VolumeList::volumeDelete(const VolumeBase* volume) {
    if (contains(volume))
        remove(volume);
}

// implementation of DataInvalidationObserver interface
void VolumeList::dataAboutToInvalidate(const DataInvalidationObservable* data) {
    notifyPendingDataInvalidation();
}

} // namespace
