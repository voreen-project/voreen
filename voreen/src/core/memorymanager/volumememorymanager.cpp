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

#include "voreen/core/memorymanager/volumememorymanager.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/gpucapabilities.h"

namespace voreen {

VolumeMemoryManager::VolumeEntity::VolumeEntity(VolumeBase* volume)
    : volume_(volume)
    , numLockedUses_(0)
{
}

const std::string VolumeMemoryManager::loggerCat_("voreen.VolumeMemoryManager");

VolumeMemoryManager::VolumeMemoryManager()
    : availableMainMemory_(0)
    , availableMainMemoryInvalid_(true)
    , availableGraphicsMemory_(0)
    , availableGraphicsMemoryInvalid_(true)
{ }

VolumeMemoryManager::~VolumeMemoryManager() {
    if (!registeredVolumes_.empty()) {
        LERROR("List of registered volumes is not empty! " << registeredVolumes_.size() << " volume(s) not deregistered.");
    }
}

void VolumeMemoryManager::registerVolume(VolumeBase* v) {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);

    if (findRegisteredVolumeEntity(v) != registeredVolumes_.end()) {
        LERROR("Cannot register volume, volume has already been registered!");
        return;
    }
    
    registeredVolumes_.push_front(v);

    updateMainMemory();
    updateGraphicsMemory();
}

void VolumeMemoryManager::deregisterVolume(VolumeBase* v) {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);
    auto it = findRegisteredVolumeEntity(v);

    if (it == registeredVolumes_.end()) {
        LERROR("Cannot deregister volume, not found in volume list!");
        return;
    }

    tgtAssert(it->numLockedUses_ == 0, "Volume Representation still has locked uses.");
    registeredVolumes_.erase(it);

    updateMainMemory();
    updateGraphicsMemory();
}

bool VolumeMemoryManager::requestMainMemory(const VolumeBase* v) {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);
    LDEBUG("Requesting RAM memory");

    v = getActualVolume(v);

    // first check: if total available main memory is not sufficient, VolumeRAM cannot be created
    if(VoreenApplication::app() && (VoreenApplication::app()->getCpuRamLimit() < getMemoryRequirement(v)))
        return false;

    size_t requiredMemory = getMemoryRequirement(v);

    // perform first check
    bool memoryCheck = (requiredMemory <= getAvailableMainMemory());

    // find the requested volume to avoid removing its representations
    auto requestedVolume = findRegisteredVolumeEntity(v);
    // start at least recently used volume to free memory
    auto currentVolume = registeredVolumes_.end() - 1;

    while (!memoryCheck) {
        LDEBUG("Not enough resources... trying to free main memory.");
        // do not remove representations from...
        if (       currentVolume == requestedVolume // the original volume
                || !(currentVolume->volume_->hasRepresentation<VolumeRAM>()) // ones that do not have a RAM representation
                || currentVolume->numLockedUses_ > 0) { // ones whose RAM representation cannot be deleted currently.
            // if every volume has been checked: break
            if (currentVolume == registeredVolumes_.begin())
                break;
            else {
                currentVolume--;
                continue;
            }
        }

        // check all representations of the current volume and see if there is one that can be converted to VolumeRAM before removing it
        if (currentVolume->volume_->canConvertToRepresentation<VolumeRAM>()) {
            currentVolume->volume_->removeRepresentation<VolumeRAM>();
            LDEBUG("Removed one RAM representation");
        }

        // update check
        memoryCheck = (requiredMemory <= getAvailableMainMemory());

        // if we have checked every volume: break
        if (currentVolume == registeredVolumes_.begin())
            break;

        // check the next volume
        currentVolume--;
    }

    LDEBUG("Main memory check result: " << memoryCheck);

    return memoryCheck;
}

bool VolumeMemoryManager::requestMainMemory(size_t requiredMemory) {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);
    LDEBUG("Requesting main memory");

    // first check: if total available main memory is not sufficient, the memory request cannot be fulfilled
    if(VoreenApplication::app() && (VoreenApplication::app()->getCpuRamLimit() < requiredMemory))
        return false;

    // perform first check
    bool memoryCheck = (requiredMemory <= getAvailableMainMemory());

    // start at least recently used volume to free memory
    auto currentVolume = registeredVolumes_.end() - 1;

    while (!memoryCheck) {
        LDEBUG("Not enough resources... trying to free main memory.");
        // do not remove representations from volumes that do not have the representation
        if (
                !(currentVolume->volume_->hasRepresentation<VolumeRAM>())
                || currentVolume->numLockedUses_ > 0) { // ... or whose RAM representation cannot be deleted currently.
            // if every volume has been checked: break
            if (currentVolume == registeredVolumes_.begin())
                break;
            else {
                currentVolume--;
                continue;
            }
        }

        // check all representations of the current volume and see if there is one that can be converted to VolumeRAM before removing it
        if (currentVolume->volume_->canConvertToRepresentation<VolumeRAM>()) {
            currentVolume->volume_->removeRepresentation<VolumeRAM>();
            LDEBUG("Removed one RAM representation");
        }

        // update check
        memoryCheck = (requiredMemory <= getAvailableMainMemory());

        // if we have checked every volume: break
        if (currentVolume == registeredVolumes_.begin())
            break;

        // check the next volume
        currentVolume--;
    }

    LDEBUG("Main memory check result: " << memoryCheck);

    return memoryCheck;
}

size_t VolumeMemoryManager::getMemoryRequirement(const VolumeBase* v) const {
    return v->getNumVoxels() * v->getBytesPerVoxel();
}

bool VolumeMemoryManager::requestGraphicsMemory(const VolumeBase* v) {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);
    LDEBUG("Requesting GPU memory");

    v = getActualVolume(v);

    // check if OpenGL types are defined, i.e., if the data type is compatible to an OpenGL texture
    if (!v->getOpenGLInternalFormat() || !v->getOpenGLFormat() || !v->getOpenGLType())
        return false;

    // first test: check if the dimensions of the volume are supported using GL_MAX_TEXTURE_SIZE
    if (tgt::max(v->getDimensions()) > static_cast<size_t>(GpuCaps.getMax3DTextureSize()))
        return false;

    // compute required memory
    size_t requiredMemory = getMemoryRequirement(v);
    // FIXME: textures seem to take up some additional space (overhead)
    requiredMemory += requiredMemory / 10;

    // check if the total GPU memory would be sufficient to upload the texture
    int totalGpuMemory = GpuCaps.retrieveTotalTextureMemory();
    if (totalGpuMemory > -1) { // check if the memory could be retrieved
        if (static_cast<int>(requiredMemory / 1024) >= totalGpuMemory) // check is performed in kilobytes
            return false;
    }

    // check if the total memory of the application settings would be sufficient
    size_t totalAvailable = VoreenApplication::app() ? VoreenApplication::app()->getGpuMemoryLimit() : 0;
    if (requiredMemory > totalAvailable)
        return false;

    // actual check: use proxy texture to check if the GPU allows the texture upload and check the available memory
    bool memoryCheck = (requiredMemory <= getAvailableGraphicsMemory()) && checkProxyTexture(v);

    // find the requested volume to avoid removing its representations
    auto requestedVolume = findRegisteredVolumeEntity(v);
    // start at least recently used volume to free memory
    auto currentVolume = registeredVolumes_.end() - 1;

    while (!memoryCheck) {
        LDEBUG("Not enough resources... trying to free GPU memory.");
        // do not remove representations from the request volumes or from volumes that do not have the representation
        if (currentVolume == requestedVolume || !(currentVolume->volume_->hasRepresentation<VolumeGL>())) {
            // if every volume has been checked: break
            if (currentVolume == registeredVolumes_.begin())
                break;
            else {
                currentVolume--;
                continue;
            }
        }

        // check all representations of the current volume and see if there is one that can be converted to VolumeGL before removing it
        if (currentVolume->volume_->canConvertToRepresentation<VolumeGL>()) {
            currentVolume->volume_->removeRepresentation<VolumeGL>();
            LDEBUG("Removed one GL representation");
        }

        // update check
        memoryCheck = (requiredMemory <= getAvailableGraphicsMemory()) && checkProxyTexture(v);

        // if we have checked every volume: break
        if (currentVolume == registeredVolumes_.begin())
            break;

        // check the next volume
        currentVolume--;
    }

    LDEBUG("GPU Memory check result: " << memoryCheck);

    return memoryCheck;
}

void VolumeMemoryManager::notifyUse(const VolumeBase* v, bool locked) {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);

    v = getActualVolume(v);

    auto it = findRegisteredVolumeEntity(v);

    if (it == registeredVolumes_.end()) {
        LERROR("Notifying use for unregistered volume!");
        return;
    }

    if (locked) {
        it->numLockedUses_++;
    }

    // remove and put to front
    VolumeEntity entity = *it;    // Backup entity.
    registeredVolumes_.erase(it); // Invalidates the iterator.
    registeredVolumes_.push_front(entity);
}

void VolumeMemoryManager::notifyLockedRelease(const VolumeBase* v) {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);

    v = getActualVolume(v);

    auto it = findRegisteredVolumeEntity(v);

    if (it == registeredVolumes_.end()) {
        LERROR("Notifying use for unregistered volume!");
        return;
    }

    tgtAssert(it->numLockedUses_ > 0, "No locked uses left");
    it->numLockedUses_--;
}

void VolumeMemoryManager::updateMainMemory() {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);
    availableMainMemoryInvalid_ = true;
}

size_t VolumeMemoryManager::getAvailableMainMemory() const {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);

    if (!availableMainMemoryInvalid_)
        return availableMainMemory_;

    // if there is no application we cannot query the available main memory
    if (!VoreenApplication::app())
        return 0;

    size_t totalMemory = VoreenApplication::app()->getCpuRamLimit();

    // iterate over volumes and add to used memory each time a VolumeRAM is found
    size_t usedMemory = 0;

    for (auto it = registeredVolumes_.begin(); it != registeredVolumes_.end(); ++it) {
        if (it->volume_->hasRepresentation<VolumeRAM>())
            usedMemory += getMemoryRequirement(it->volume_);
    }

    if (usedMemory > totalMemory)   // may occur if the application settings have changed
        return 0;
    else
        return totalMemory - usedMemory;
}

void VolumeMemoryManager::updateGraphicsMemory() {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);
    availableGraphicsMemoryInvalid_ = true;
}

size_t VolumeMemoryManager::getAvailableGraphicsMemory() const {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);

    if (!availableGraphicsMemoryInvalid_)
        return availableGraphicsMemory_;

    // if there is no application we cannot query the available main memory
    if (!VoreenApplication::app())
        return 0;

    size_t totalMemory = VoreenApplication::app()->getGpuMemoryLimit();

    // iterate over volumes and add to used memory each time a VolumeGL is found
    size_t usedMemory = 0;

    for (auto it = registeredVolumes_.begin(); it != registeredVolumes_.end(); ++it) {
        if (it->volume_->hasRepresentation<VolumeGL>()) {
            size_t requiredMemory = getMemoryRequirement(it->volume_);
            // FIXME: textures seem to take up some additional space (overhead)
            requiredMemory += requiredMemory / 10;
            usedMemory += requiredMemory;
        }
    }

    if (usedMemory > totalMemory)   // may occur if the application settings have changed
        return 0;
    else
        return totalMemory - usedMemory;
}

bool VolumeMemoryManager::checkProxyTexture(const VolumeBase* v) {
    boost::lock_guard<boost::recursive_mutex> lock(vmmMutex_);

    tgt::svec3 volDim = v->getDimensions();

    // call glTexImage3D to test the OpenGL implementation
    glTexImage3D(GL_PROXY_TEXTURE_3D, 0, v->getOpenGLInternalFormat(), volDim.x, volDim.y, volDim.z, 0, v->getOpenGLFormat(), v->getOpenGLType(), 0);

    // check if the call was successful
    GLint width;
    glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0, GL_TEXTURE_WIDTH, &width);
    if (width == 0)
        return false;

    if (width != static_cast<int>(volDim.x)) {
        LWARNING("Width of proxy texture = " << width << ", width of volume = " << volDim.x);
        return false;
    }

    return true;
}

boost::recursive_mutex* VolumeMemoryManager::getMutex() {
    return &vmmMutex_;
}

const VolumeBase* VolumeMemoryManager::getActualVolume(const VolumeBase* v) {
    while (const VolumeDecoratorIdentity* dec = dynamic_cast<const VolumeDecoratorIdentity*>(v)) {
        v = dec->getDecorated();
    }

    return dynamic_cast<const Volume*>(v);
}

std::deque<VolumeMemoryManager::VolumeEntity>::iterator VolumeMemoryManager::findRegisteredVolumeEntity(const VolumeBase* v) {
    return std::find_if(registeredVolumes_.begin(), registeredVolumes_.end(), [v] (VolumeEntity& ref) {
            return ref.volume_ == v;
            });
}

VolumeRAMRepresentationLock::VolumeRAMRepresentationLock(const VolumeBase* volume)
    : volume_(volume)
    , representation_(nullptr)
{
    tgtAssert(volume_, "volume was null");

    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().notifyUse(volume_, true);
    representation_ = volume_->getRepresentation<VolumeRAM>();
}
VolumeRAMRepresentationLock::VolumeRAMRepresentationLock(const VolumeRAMRepresentationLock& other)
    : VolumeRAMRepresentationLock(other.volume_)
{
}
VolumeRAMRepresentationLock::~VolumeRAMRepresentationLock() {
    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().notifyLockedRelease(volume_);
}

VolumeRAMRepresentationLock& VolumeRAMRepresentationLock::operator=(const VolumeBase* volume) {
    tgtAssert(volume, "volume was null");

    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().notifyLockedRelease(volume_);

    volume_ = volume;

    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().notifyUse(volume_, true);

    representation_ = volume->getRepresentation<VolumeRAM>();

    return *this;
}
VolumeRAMRepresentationLock& VolumeRAMRepresentationLock::operator=(const VolumeRAMRepresentationLock& other) {
    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().notifyLockedRelease(volume_);

    volume_ = other.volume_;

    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().notifyUse(volume_, true);

    representation_ = other.representation_;

    return *this;
}

const VolumeRAM* VolumeRAMRepresentationLock::operator->() const {
    return representation_;
}
const VolumeRAM* VolumeRAMRepresentationLock::operator*() const {
    return representation_;
}

} // namespace voreen
