/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "voreen/core/datastructures/volume/volume.h"

// hashing
#include "voreen/core/datastructures/volume/volumehash.h"

// representations
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumegl.h"
#include "voreen/core/datastructures/volume/volumedisk.h"

// application (for progress bar when loading a volume)
#include "voreen/core/voreenapplication.h"
#include "voreen/core/io/progressbar.h"

// serialization
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumeserializer.h"

namespace voreen {

const std::string Volume::loggerCat_("voreen.Volume");

//---- constructor and destructor ----

Volume::Volume(VolumeRepresentation* const volume, const tgt::vec3& spacing, const tgt::vec3& offset, const tgt::mat4& transformation)
{
    format_ = volume->getFormat();
    baseType_ = volume->getBaseType();
    // const cast is safe because we get a non-const volume originally
    determineOpenGLTypes(format_);

    setSpacing(spacing);
    setOffset(offset);
    setPhysicalToWorldMatrix(transformation);

    addRepresentation(volume);

    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().registerVolume(this);
}

Volume::Volume(VolumeRepresentation* const volume, const VolumeBase* vh)
{
    format_ = volume->getFormat();
    baseType_ = volume->getBaseType();
    determineOpenGLTypes(format_);

    std::vector<std::string> keys = vh->getMetaDataKeys();
    for(size_t i=0; i<keys.size(); i++) {
        const MetaDataBase* md = vh->getMetaData(keys[i]);
        if(md) {
            metaData_.addMetaData(keys[i], md->clone());
        }
    }

    addRepresentation(volume);

    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().registerVolume(this);
}

Volume::Volume(VolumeRepresentation* const volume, const MetaDataContainer* mdc) {

    format_ = volume->getFormat();
    baseType_ = volume->getBaseType();
    determineOpenGLTypes(format_);

    std::vector<std::string> keys = mdc->getKeys();
    for(size_t i=0; i<keys.size(); i++) {
        const MetaDataBase* md = mdc->getMetaData(keys[i]);
        if(md) {
            metaData_.addMetaData(keys[i], md->clone());
        }
    }

    addRepresentation(volume);

    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().registerVolume(this);
}

Volume::Volume(VolumeRepresentation* const volume, const MetaDataContainer* mdc, const std::set<VolumeDerivedData*>& derivedData) {

    format_ = volume->getFormat();
    baseType_ = volume->getBaseType();
    determineOpenGLTypes(format_);

    std::vector<std::string> keys = mdc->getKeys();
    for(size_t i=0; i<keys.size(); i++) {
        const MetaDataBase* md = mdc->getMetaData(keys[i]);
        if(md) {
            metaData_.addMetaData(keys[i], md->clone());
        }
    }

    for(std::set<VolumeDerivedData*>::const_iterator i = derivedData.begin(); i != derivedData.end(); i++)
       addDerivedData(*i);

    addRepresentation(volume);

    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().registerVolume(this);
}

Volume::~Volume() {
    notifyDelete();
    stopRunningThreads();
    deleteAllRepresentations(false);

    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().deregisterVolume(this);
}


//---- hash ----

void Volume::setHash(const std::string& hash) const {
    addDerivedDataInternal<VolumeHash>(new VolumeHash(hash));
}


//---- meta data ----

const MetaDataContainer& Volume::getMetaDataContainer() const {
    return metaData_;
}

MetaDataContainer& Volume::getMetaDataContainer() {
    return metaData_;
}


std::vector<std::string> Volume::getMetaDataKeys() const {
    return metaData_.getKeys();
}

const MetaDataBase* Volume::getMetaData(const std::string& key) const {
    return metaData_.getMetaData(key);
}

bool Volume::hasMetaData(const std::string& key) const {
    return metaData_.hasMetaData(key);
}

void Volume::setModality(const Modality& modality) {
    setMetaDataValue<StringMetaData>(META_DATA_NAME_MODALITY, modality.getName());
}

void Volume::setRealWorldMapping(const RealWorldMapping& rwm) {
    setMetaDataValue<RealWorldMappingMetaData>(META_DATA_NAME_REAL_WORLD_MAPPING, rwm);
}

void Volume::setTimestep(float timestep) {
    setMetaDataValue<FloatMetaData>(META_DATA_NAME_TIMESTEP, timestep);
}

void Volume::setSpacing(const tgt::vec3& spacing) {
    setMetaDataValue<Vec3MetaData>(META_DATA_NAME_SPACING, spacing);
}

void Volume::setOffset(const tgt::vec3& offset) {
    setMetaDataValue<Vec3MetaData>(META_DATA_NAME_OFFSET, offset);
}

void Volume::setPhysicalToWorldMatrix(const tgt::mat4& transformationMatrix) {
    return setMetaDataValue<Mat4MetaData>(META_DATA_NAME_TRANSFORMATION, transformationMatrix);
}


//---- representation data ----

size_t Volume::getNumRepresentations() const {
    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    return representations_.size();
}

const VolumeRepresentation* Volume::getRepresentation(size_t i) const {
    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);

    // do not notify memory manager that this volume is in use, since this is an internal method!
    //if (VolumeMemoryManager::isInited())
    //    VolumeMemoryManager::getRef().notifyUse(const_cast<VolumeBase*>(this));

    return representations_[i];
}

const VolumeRepresentation* Volume::useConverter(const RepresentationConverterBase* converter) const {

    // volumes must lock the memory manager if they use it to prevent deadlocks in multi-threading
    boost::recursive_mutex* vmmMutex = 0;
    if (VolumeMemoryManager::isInited())
        vmmMutex = VolumeMemoryManager::getRef().getMutex();
    VolumeLockGuard vmmGuard(vmmMutex);

    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    for(size_t i=0; i < representations_.size(); i++) {
        if(converter->canConvert(getRepresentation(i))) {
            VolumeRepresentation* rep = converter->convert(getRepresentation(i));

            if(rep) {

                if (VolumeMemoryManager::isInited()) {
                    if (dynamic_cast<VolumeRAM*>(rep))
                        VolumeMemoryManager::getRef().updateMainMemory();
                    else if (dynamic_cast<VolumeGL*>(rep))
                        VolumeMemoryManager::getRef().updateGraphicsMemory();
                }

                representations_.push_back(rep);
                return rep;
            }
        }
    }
    return 0;
}

void Volume::addRepresentation(VolumeRepresentation* rep) {

    // volumes must lock the memory manager if they use it to prevent deadlocks in multi-threading
    boost::recursive_mutex* vmmMutex = 0;
    if (VolumeMemoryManager::isInited())
        vmmMutex = VolumeMemoryManager::getRef().getMutex();
    VolumeLockGuard vmmGuard(vmmMutex);

    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    //TODO: check for duplicates using RTI

    // check if the memory manager has to be notified to update its main memory or graphics memory (is done lazy, so we can call it before pushing back the representation)
    if (VolumeMemoryManager::isInited()) {
        if (dynamic_cast<VolumeRAM*>(rep))
            VolumeMemoryManager::getRef().updateMainMemory();
        else if (dynamic_cast<VolumeGL*>(rep))
            VolumeMemoryManager::getRef().updateGraphicsMemory();
    }

    representations_.push_back(rep);
}

void Volume::removeRepresentation(size_t i) {

    // volumes must lock the memory manager if they use it to prevent deadlocks in multi-threading
    boost::recursive_mutex* vmmMutex = 0;
    if (VolumeMemoryManager::isInited())
        vmmMutex = VolumeMemoryManager::getRef().getMutex();
    VolumeLockGuard vmmGuard(vmmMutex);

    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    if (i >= representations_.size())
        return;

    notifyRepresentationDelete(representations_.at(i));

    // we can safely notify the memory manager here as it updates the main memory and graphics memory lazy!
    if (VolumeMemoryManager::isInited()) {
        if (dynamic_cast<VolumeRAM*>(representations_.at(i)))
            VolumeMemoryManager::getRef().updateMainMemory();
        else if (dynamic_cast<VolumeGL*>(representations_.at(i)))
            VolumeMemoryManager::getRef().updateGraphicsMemory();
    }

    stopRunningThreads();
    delete representations_.at(i);
    representations_.erase(representations_.begin() + i);
}

void Volume::deleteAllRepresentations(bool notifyObserver /*= true*/) {

    // volumes must lock the memory manager if they use it to prevent deadlocks in multi-threading
    boost::recursive_mutex* vmmMutex = 0;
    if (VolumeMemoryManager::isInited())
        vmmMutex = VolumeMemoryManager::getRef().getMutex();
    VolumeLockGuard vmmGuard(vmmMutex);

    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    stopRunningThreads();

    bool hasVolumeRAM = VolumeBase::hasRepresentation<VolumeRAM>();
    bool hasVolumeGL = VolumeBase::hasRepresentation<VolumeGL>();

    // notify observers we will delete the representations
    if (notifyObserver) {
        for(size_t i=0; i < representations_.size(); i++) {
            notifyRepresentationDelete(representations_[i]);
        }
    }

    while(!representations_.empty()) {
        delete representations_.back();
        representations_.pop_back();
    }

    if (VolumeMemoryManager::isInited()) {
        if (hasVolumeRAM)
            VolumeMemoryManager::getRef().updateMainMemory();

        if (hasVolumeGL)
            VolumeMemoryManager::getRef().updateGraphicsMemory();
    }
}

void Volume::releaseAllRepresentations() {

    // volumes must lock the memory manager if they use it to prevent deadlocks in multi-threading
    boost::recursive_mutex* vmmMutex = 0;
    if (VolumeMemoryManager::isInited())
        vmmMutex = VolumeMemoryManager::getRef().getMutex();
    VolumeLockGuard vmmGuard(vmmMutex);

    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);

    bool hasVolumeRAM = VolumeBase::hasRepresentation<VolumeRAM>();
    bool hasVolumeGL = VolumeBase::hasRepresentation<VolumeGL>();

    // notify observers as if we have deleted the representations
    for(size_t i=0; i < representations_.size(); i++) {
        notifyRepresentationDelete(representations_[i]);
    }

    stopRunningThreads();
    representations_.clear();

    if (VolumeMemoryManager::isInited()) {
        if (hasVolumeRAM)
            VolumeMemoryManager::getRef().updateMainMemory();

        if (hasVolumeGL)
            VolumeMemoryManager::getRef().updateGraphicsMemory();
    }
}


void Volume::setVolumeData(VolumeRAM* const volume) {
    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    if(!volume) {
        LERROR("Tried to set null volume!");
        tgtAssert(false, "Tried to set null volume!");
        return;
    }

    if(VolumeBase::hasRepresentation<VolumeRAM>()) {
        const VolumeRAM* v = getRepresentation<VolumeRAM>();

        if(v != volume) {
            if(v->getDimensions() != volume->getDimensions()) {
                LERROR("Tried to set volume with different dimensions!");
                tgtAssert(false, "Tried to set volume with different dimensions!");
                return;
            }

            deleteAllRepresentations();
            addRepresentation(volume);
        }
    }
    else {
        addRepresentation(volume);
        makeRepresentationExclusive<VolumeRAM>();
    }
}

bool Volume::reloadVolume() {
    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);

    ProgressBar* progressDialog = VoreenApplication::app()->createProgressDialog();
    if (progressDialog) {
        progressDialog->setTitle("Loading volume");
        progressDialog->setProgressMessage("Loading volume ...");
    }
    VolumeSerializerPopulator populator(progressDialog);

    // try to load volume from origin
    VolumeBase* handle = 0;
    try {
        handle = populator.getVolumeSerializer()->read(origin_);
    }
    catch (tgt::FileException& e) {
        LWARNING(e.what());
    }
    catch (std::bad_alloc&) {
        LWARNING("std::Error BAD ALLOCATION");
    }
    delete progressDialog;
    progressDialog = 0;

    if (!handle || !handle->getRepresentation<VolumeRAM>()) {
        delete handle;
        return false;
    }

    Volume* vh = static_cast<Volume*>(handle);
    if(VolumeBase::hasRepresentation<VolumeRAM>()) {
        deleteAllRepresentations();

        addRepresentation(vh->getWritableRepresentation<VolumeRAM>());
        vh->releaseAllRepresentations();
        delete handle;
    }

    // inform observers
    notifyChanged();
    return true;
}

void Volume::invalidate() {
    stopRunningThreads();
    clearFinishedThreads();
    for(std::set<VolumeDerivedData*>::iterator it = derivedData_.begin(); it != derivedData_.end(); it++)
        delete (*it);
    derivedData_.clear();

    // volumes must lock the memory manager if they use it to prevent deadlocks in multi-threading
    boost::recursive_mutex* vmmMutex = 0;
    if (VolumeMemoryManager::isInited())
        vmmMutex = VolumeMemoryManager::getRef().getMutex();
    VolumeLockGuard vmmGuard(vmmMutex);

    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);

    if(VolumeBase::hasRepresentation<VolumeGL>()) {
        VolumeBase::removeRepresentation<VolumeGL>();
        if (VolumeMemoryManager::isInited())
            VolumeMemoryManager::getRef().updateGraphicsMemory();
    }
    if(VolumeBase::hasRepresentation<VolumeRAM>()) {
        getRepresentation<VolumeRAM>()->invalidate();
        if (VolumeMemoryManager::isInited())
            VolumeMemoryManager::getRef().updateMainMemory();
    }
    if(VolumeBase::hasRepresentation<VolumeDisk>()) {
        getRepresentation<VolumeDisk>()->invalidate();
    }
}

void Volume::notifyRepresentationDelete(const VolumeRepresentation* rep) {
    std::vector<VolumeObserver*> observers = Observable<VolumeObserver>::getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->volumeRepresentationDelete(this, rep);
}

// ---- basic volume information ----

void Volume::determineOpenGLTypes(std::string format) {
    if (format_ == "uint8") {
        glFormat_= GL_RED;
        glInternalFormat_= GL_R8;
        glType_ = GL_UNSIGNED_BYTE;
    }
    else if (format_ == "int8") {
        glFormat_= GL_RED;
        glInternalFormat_= GL_R8_SNORM;
        glType_ = GL_BYTE;
    }
    else if (format_ == "uint16") {
        glFormat_= GL_RED;
        glInternalFormat_= GL_R16;
        glType_ = GL_UNSIGNED_SHORT;
    }
    else if (format_ == "int16") {
        glFormat_= GL_RED;
        glInternalFormat_= GL_R16_SNORM;
        glType_ = GL_SHORT;
    }
    else if (format_ == "uint32") {
        glFormat_= GL_RED;
        glInternalFormat_= GL_R32F;
        glType_ = GL_UNSIGNED_INT;
    }
    else if (format_ == "int32") {
        glFormat_= GL_RED;
        glInternalFormat_= GL_R32F;
        glType_ = GL_INT;
    }
    else if (format_ == "float") {
        glFormat_= GL_RED;
        glInternalFormat_= GL_R32F;
        glType_ = GL_FLOAT;
    }
    else if (format_ == "Vector2(uint8)") {
        glFormat_= GL_RG;
        glInternalFormat_= GL_RG8;
        glType_ = GL_UNSIGNED_BYTE;
    }
    else if (format_ == "Vector2(int8)") {
        glFormat_= GL_RG;
        glInternalFormat_= GL_RG8_SNORM;
        glType_ = GL_BYTE;
    }
    else if (format_ == "Vector2(uint16)") {
        glFormat_= GL_RG;
        glInternalFormat_= GL_RG16;
        glType_ = GL_UNSIGNED_SHORT;
    }
    else if (format_ == "Vector2(int16)") {
        glFormat_= GL_RG;
        glInternalFormat_= GL_RG16_SNORM;
        glType_ = GL_SHORT;
    }
    else if (format_ == "Vector2(uint32)") {
        glFormat_= GL_RG;
        glInternalFormat_= GL_RG32F;
        glType_ = GL_UNSIGNED_INT;
    }
    else if (format_ == "Vector2(int32)") {
        glFormat_= GL_RG;
        glInternalFormat_= GL_RG32F;
        glType_ = GL_INT;
    }
    else if (format_ == "Vector2(float)") {
        glFormat_= GL_RG;
        glInternalFormat_= GL_RG32F;
        glType_ = GL_FLOAT;
    }
    else if (format_ == "Vector3(uint8)") {
        glFormat_= GL_RGB;
        glInternalFormat_= GL_RGB8;
        glType_ = GL_UNSIGNED_BYTE;
    }
    else if (format_ == "Vector3(int8)") {
        glFormat_= GL_RGB;
        glInternalFormat_= GL_RGB8_SNORM;
        glType_ = GL_BYTE;
    }
    else if (format_ == "Vector3(uint16)") {
        glFormat_= GL_RGB;
        glInternalFormat_= GL_RGB16;
        glType_ = GL_UNSIGNED_SHORT;
    }
    else if (format_ == "Vector3(int16)") {
        glFormat_= GL_RGB;
        glInternalFormat_= GL_RGB16_SNORM;
        glType_ = GL_SHORT;
    }
    else if (format_ == "Vector3(uint32)") {
        glFormat_= GL_RGB;
        glInternalFormat_= GL_RGB;
        glType_ = GL_UNSIGNED_INT;
    }
    else if (format_ == "Vector3(int32)") {
        glFormat_= GL_RGB;
        glInternalFormat_= GL_RGB32F;
        glType_ = GL_INT;
    }
    else if (format_ == "Vector3(float)") {
        glFormat_= GL_RGB;
        glInternalFormat_= GL_RGB32F;
        glType_ = GL_FLOAT;
    }
    else if (format_ == "Vector4(uint8)") {
        glFormat_= GL_RGBA;
        glInternalFormat_= GL_RGBA8;
        glType_ = GL_UNSIGNED_BYTE;
    }
    else if (format_ == "Vector4(int8)") {
        glFormat_= GL_RGBA;
        glInternalFormat_= GL_RGBA8_SNORM;
        glType_ = GL_BYTE;
    }
    else if (format_ == "Vector4(uint16)") {
        glFormat_= GL_RGBA;
        glInternalFormat_= GL_RGBA16;
        glType_ = GL_UNSIGNED_SHORT;
    }
    else if (format_ == "Vector4(int16)") {
        glFormat_= GL_RGBA;
        glInternalFormat_= GL_RGBA16_SNORM;
        glType_ = GL_SHORT;
    }
    else if (format_ == "Vector4(uint32)") {
        glFormat_= GL_RGBA;
        glInternalFormat_= GL_RGBA;
        glType_ = GL_UNSIGNED_INT;
    }
    else if (format_ == "Vector4(int32)") {
        glFormat_= GL_RGBA;
        glInternalFormat_= GL_RGBA32F;
        glType_ = GL_INT;
    }
    else if (format_ == "Vector4(float)") {
        glFormat_= GL_RGBA;
        glInternalFormat_= GL_RGBA32F;
        glType_ = GL_FLOAT;
    }
    else {
        // unknown or not supported format (e.g., double) -> set values to 0 to produce OpenGL error
        glInternalFormat_ = 0;
        glFormat_ = 0;
        glType_ = 0;
        LDEBUG("Could not determine OpenGL types!");
    }
}

//---------------------------------------------------------------

void oldVolumePosition(Volume* vh) {
    //correct old spacing:
    tgt::vec3 sp = vh->getSpacing();
    tgt::vec3 cubeSize = sp * tgt::vec3(vh->getDimensions());

    float scale = 2.0f / max(cubeSize);
    sp *= scale;
    vh->setSpacing(sp);

    //set origin to center volume:
    cubeSize = sp * tgt::vec3(vh->getDimensions());
    vh->setOffset(-cubeSize/2.0f);
}

void centerVolume(Volume* vh) {
    tgt::vec3 sp = vh->getSpacing();
    tgt::vec3 cubeSize = sp * tgt::vec3(vh->getDimensions());

    //set origin to center volume:
    vh->setOffset(-cubeSize/2.0f);
}

} // namespace
