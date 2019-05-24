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

#include "voreen/core/datastructures/volume/volumedecorator.h"

#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/volume/volumehash.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/utils/hashing.h"

using std::string;
using tgt::vec3;
using tgt::mat4;

namespace voreen {

//---- constructor / destructor ----

VolumeDecoratorIdentity::VolumeDecoratorIdentity(const VolumeBase* vhb)
    : VolumeBase()
    , VolumeObserver()
    , base_(vhb)
{
    tgtAssert(vhb, "null pointer passed as VolumeBase");
    vhb->Observable<VolumeObserver>::addObserver(this);
}

VolumeDecoratorIdentity::~VolumeDecoratorIdentity() {
    notifyDelete();
    stopRunningThreads();
}

const VolumeBase* VolumeDecoratorIdentity::getDecorated() const {
    return base_;
}


//---- basic volume information ----
std::string VolumeDecoratorIdentity::getFormat() const {
    if (base_)
        return base_->getFormat();
    else
        return format_;
}

std::string VolumeDecoratorIdentity::getBaseType() const {
    if (base_)
        return base_->getBaseType();
    else
        return baseType_;
}

GLint VolumeDecoratorIdentity::getOpenGLInternalFormat() const {
    if (base_)
        return base_->getOpenGLInternalFormat();
    else
        return glInternalFormat_;
}

GLenum VolumeDecoratorIdentity::getOpenGLFormat() const {
    if (base_)
        return base_->getOpenGLFormat();
    else
        return glFormat_;
}

GLenum VolumeDecoratorIdentity::getOpenGLType() const {
    if (base_)
        return base_->getOpenGLType();
    else
        return glType_;
}


//---- representations ----

size_t VolumeDecoratorIdentity::getNumRepresentations() const {
    if (base_)
        return base_->getNumRepresentations();
    else
        return 0;
}

const VolumeRepresentation* VolumeDecoratorIdentity::getRepresentation(size_t i) const {
    if (base_)
        return base_->getRepresentation(i);
    else
        return 0;
}

void VolumeDecoratorIdentity::removeRepresentation(size_t i) {
    stopRunningThreads();

    if (base_)
        const_cast<VolumeBase*>(base_)->removeRepresentation(i);
    else
        LERROR("Trying to remove representation from VolumeDecorator without decorated volume!");
}

const VolumeRepresentation* VolumeDecoratorIdentity::useConverter(const RepresentationConverterBase* converter) const {
    if (base_)
        return base_->useConverter(converter);
    else
        return 0;
}

void VolumeDecoratorIdentity::addRepresentation(VolumeRepresentation* rep) {
    if (base_)
        const_cast<VolumeBase*>(base_)->addRepresentation(rep);
    else
        LERROR("Trying to add representation to VolumeDecorator without decorated volume!");
}


//---- meta data -----

std::vector<std::string> VolumeDecoratorIdentity::getMetaDataKeys() const {
    if (base_)
        return base_->getMetaDataKeys();
    else {
        std::vector<std::string> emptyVector;
        return emptyVector;
    }
}

const MetaDataBase* VolumeDecoratorIdentity::getMetaData(const std::string& key) const {
    if (base_)
        return base_->getMetaData(key);
    else
        return 0;
}

bool VolumeDecoratorIdentity::hasMetaData(const std::string& key) const {
    if (base_)
        return base_->hasMetaData(key);
    else
        return false;
}

Modality VolumeDecoratorIdentity::getModality() const {
    if (base_)
        return base_->getModality();
    else
        return Modality("unknown");
}


// ---- volume observer ----

void VolumeDecoratorIdentity::volumeDelete(const VolumeBase* source) {
    // note: this is also called in volumeDataDelete (-> change volumeDataDelete when this implementation changes, if necessary)
    tgtAssert(source == base_, "Notified by other volume than decorated!");

    // stop derived data threads
    stopRunningThreads();
    clearFinishedThreads();
    clearDerivedData();

    // not ideal, but the best we can do here: notify observers about deletion of data
    base_->Observable<VolumeObserver>::removeObserver(this);
    std::vector<VolumeObserver*> observers = Observable<VolumeObserver>::getObservers();
    for (size_t i = 0; i < observers.size(); ++i)
        observers[i]->volumeDataDelete(this);
    base_ = 0;
}

void VolumeDecoratorIdentity::volumeChange(const VolumeBase* /*source*/) {
    notifyChanged();
    // stop derived data threads
    stopRunningThreads();
    clearFinishedThreads();
    clearDerivedData();
}

void VolumeDecoratorIdentity::volumeRepresentationDelete(const VolumeBase* volume, const VolumeRepresentation* rep) {
    // if a representation of the original volume has been removed, we should notify our observers, too
    tgtAssert(volume == base_, "Notified not by decorated volume");

    std::vector<VolumeObserver*> observers = Observable<VolumeObserver>::getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->volumeRepresentationDelete(this, rep);
}

void VolumeDecoratorIdentity::volumeDataDelete(const VolumeBase* source) {
    // we have to stop threads, notify observers, and set base_ to 0 -> we just call volumeDelete, which does exactly those things (atm)
    volumeDelete(source);
}

//-------------------------------------------------------------------------------------------------

VolumeDecoratorReplace::VolumeDecoratorReplace(const VolumeBase* vhb,
    const std::string& key, MetaDataBase* value, bool keepDerivedData)
    : VolumeDecoratorIdentity(vhb)
    , key_(key)
    , value_(value)
{
    tgtAssert(key != "", "empty key passed");
    tgtAssert(value, "null pointer passed as value");

    // create volume hash by concatenating hash of underlying volume with the replace item
    if (vhb->hasDerivedData<VolumeHash>()) {
        VolumeHash* newHash = new VolumeHash();
        newHash->setHash(VoreenHash::getHash(vhb->getHash() + "-" + value->toString()));
        addDerivedData(newHash);
    }

    // copy over derived data (TODO: handle other derived data types)
    if (keepDerivedData) {
        if (vhb->hasDerivedData<VolumeHistogramIntensity>())
            addDerivedData(new VolumeHistogramIntensity(*vhb->getDerivedData<VolumeHistogramIntensity>()));

        if (vhb->hasDerivedData<VolumeMinMax>())
            addDerivedData(new VolumeMinMax(*vhb->getDerivedData<VolumeMinMax>()));

        if (vhb->hasDerivedData<VolumeMinMaxMagnitude>())
            addDerivedData(new VolumeMinMaxMagnitude(*vhb->getDerivedData<VolumeMinMaxMagnitude>()));
    }

    // copy over volume origin
    setOrigin(vhb->getOrigin());
}

std::vector<std::string> VolumeDecoratorReplace::getMetaDataKeys() const {
    if (base_) {
        std::vector<std::string> keys = base_->getMetaDataKeys();

        if(!base_->hasMetaData(key_))
            keys.push_back(key_);

        return keys;
    }
    else {
        std::vector<std::string> emptyVector;
        return emptyVector;
    }
}

const MetaDataBase* VolumeDecoratorReplace::getMetaData(const std::string& key) const {
    if (!base_)
        return 0;

    if(key == key_)
        return value_;
    else
        return base_->getMetaData(key);
}

bool VolumeDecoratorReplace::hasMetaData(const std::string& key) const {
    if (!base_)
        return false;

    if(key == key_)
        return true;
    else
        return base_->hasMetaData(key);
}

MetaDataBase* VolumeDecoratorReplace::getValue() const {
    return value_;
}

void VolumeDecoratorReplace::setValue(MetaDataBase* value) {
    notifyPendingDataInvalidation();
    delete value_;
    value_ = value;
}

} // namespace
