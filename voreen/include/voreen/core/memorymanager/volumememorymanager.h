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

#ifndef VRN_VOLUMEMEMORYMANAGER_H
#define VRN_VOLUMEMEMORYMANAGER_H

#include "tgt/singleton.h"

#include "voreen/core/voreencoreapi.h"

#include <boost/thread.hpp>
#include <deque>

namespace voreen {

class VolumeBase;
class VolumeRAM;
class VolumeRAMRepresentationLock;

struct VolumeRef {
    VolumeBase* volume_;
    size_t numLockedUses_;

    VolumeRef(VolumeBase* volume);
    ~VolumeRef();
};

/**
 * This class provides basic memory management for volume data and its representations.
 * Each volume (i.e. VolumeBase object) has to be registered (which is managed transparently by the VolumeBase class constructor / destructor).
 * Once a representation of a Volume is requested, the memory manager checks the available memory and, if it is not sufficient, removes representations from other volumes using an LRU strategy.
 * The main memory and GPU memory / representations are thereby handled individually.
 */
class VRN_CORE_API VolumeMemoryManager : public tgt::Singleton<VolumeMemoryManager> {

public:

    /**
     *   Init VolumeMemoryManager.
     */
    VolumeMemoryManager();
    virtual ~VolumeMemoryManager();

protected:

    // VolumeBase and Volume have to call protected functions in MemoryManager and are therefore friend classes
    friend class VolumeBase;
    friend class Volume;

    // the application has to notify the memory manager if the memory settings have changed
    friend class VoreenApplication;

    // the converter from disk to GL internally converts to a VolumeRAM and should therefore use the memory manager
    friend class RepresentationConverterLoadFromDiskToGL;

    // Needs to call notify*-methods
    friend class VolumeRAMRepresentationLock;

    /**
     * Register volume for memory management (called by VolumeBase constructor).
     */
    void registerVolume(VolumeBase* v);

    /**
     * Deregister volume (called by VolumeBase destructor).
     */
    void deregisterVolume(VolumeBase* v);

    /**
     * Notifies the memory manager that the volume is used (i.e., is called by VolumeBase class if a representation is requested).
     * This will place the volume at the front of the LRU list (as it is the most recently used volume).
     *
     * @param locked Specify whether or not the RAM representation of volume should also be considered locked to main memory
     * and thus can be freed if other sources request main memory AND if there are no other sources that lock this
     * RAM representation.
     */
    void notifyUse(const VolumeBase* v, bool locked = false);

    /**
     * Notify the memory manager that the RAM representation of the specified volume is no longer locked into memory and
     * thus can be freed if other sources request main memory AND if there are no other sources that lock this RAM representation.
     */
    void notifyLockedRelease(const VolumeBase* v);

    /**
     * Request main memory, e.g. for a VolumeRAM representation (called by getRepresentation).
     *
     * @return true if the memory can be provided, false if not
     */
    bool requestMainMemory(const VolumeBase* v);

    /**
     * Request main memory, e.g. for internally converting to a VolumeRAM representation.
     *
     * Used by the Disk-to-GL converter, which internally creates a temporary VolumeRAM.
     */
    bool requestMainMemory(size_t requiredMemory);

    /**
     * Request graphics (texture) memory, e.g. for a VolumeGL representation (called by getRepresentation).
     *
     * @return true if the memory could be provided, false if not
     */
    bool requestGraphicsMemory(const VolumeBase* v);

    /**
     * Notifies the VolumeMemoryManager that the current available main memory has to be updated (e.g., when removing or adding a VolumeRAM representation).
     * Only sets an invalid flag and updates the actual memory in a lazy manner, so it can be called in several spots in the volume class!
     */
    void updateMainMemory();

    /**
     * Computes the main memory currently available for VolumeRAM representations (in bytes).
     */
    size_t getAvailableMainMemory() const;

    /**
     * Notifies the VolumeMemoryManager that the currently available GPU texture memory has to be udated (e.g., when removing or adding a VolumeGL representation).
     * Only sets the invalid flag and updates the actual memory in a lazy manner, so it can be called in several spots in the volume class!
     */
    void updateGraphicsMemory();

    /**
     * Computes the graphics memory currently available for VolumeGL representations (in bytes).
     */
    size_t getAvailableGraphicsMemory() const;

    /**
     * Uses OpenGL proxy textures to check if a volume texture could be uploaded.
     */
    bool checkProxyTexture(const VolumeBase* v);

    /// Returns the memory required for this volume (in bytes)
    size_t getMemoryRequirement(const VolumeBase* v) const;

    /// internal method for finding a decorated volume in a hierarchy of VolumeDecorators
    const VolumeBase* getActualVolume(const VolumeBase* v);

    /**
     * Returns a pointer to the mutex of the memory manager, which is used in Volume and VolumeBase class.
     *
     * Should only be used when absolutely necessary!!!
     */
    boost::recursive_mutex* getMutex();

    static const std::string loggerCat_;

    std::deque<VolumeRef> registeredVolumes_;   ///< the list of registered volumes in LRU order (front: most recently used, back: lest recently used)

    std::deque<VolumeRef>::iterator findRegisteredVolume(const VolumeBase*);

    size_t availableMainMemory_;            ///< currently available main memory for volumes (in bytes)
    bool availableMainMemoryInvalid_;       ///< does the available main memory have to be recomputed?

    size_t availableGraphicsMemory_;        ///< currently available graphics memory for volumes (in bytes)
    bool availableGraphicsMemoryInvalid_;   ///< does the available graphics memory has to be recomputed?

    mutable boost::recursive_mutex vmmMutex_;     ///< mutex for the volume memory manager
};

// Use this class to get a VolumeRAM Representation from a volume and force it to remain in main memory until
// this object is dropped.
class VolumeRAMRepresentationLock {
public:
    const VolumeRAM* operator->() const;
private:
    const VolumeBase* vol_;
    const VolumeRAM* ram_;

    VolumeRAMRepresentationLock(const VolumeBase* vol);
    ~VolumeRAMRepresentationLock();
};

} // namespace voreen


#endif //VRN_VOLUMEMEMORYMANAGER_H
