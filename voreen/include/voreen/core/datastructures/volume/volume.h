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

#ifndef VRN_VOLUME_H
#define VRN_VOLUME_H

#include "volumebase.h"

namespace voreen {

/**
 * Class for handling the actual volume data and its representations, including the volume's meta data.
 *
 * This is designed to be the only class that handles the actual data representations (although the VolumeBase class provides a basic interface which should be sufficient for data access).
 */
class VRN_CORE_API Volume : public VolumeBase {
public:

    /**
     * Constructor.
     *
     * @note No hardware specific volume data like VolumeGL are created initially. If you want
     *  to use hardware specific volume data / textures, call getRepresentation<T>() with the desired type.
     *
     * @param   volume         The volume data for this Volume.
     * @param   spacing        The distance between two neighboring voxels in all 3 dimensions (physical coordinates).
     * @param   offset         The offset of the (center of the) voxel with coordinates (0,0,0) from the origin (physical coordinates).
     * @param   transformation Is used to transform physical coordinates (defined by origin and spacing) to world coordinates.
     */
    Volume(VolumeRepresentation* const volume, const tgt::vec3& spacing, const tgt::vec3& offset, const tgt::mat4& transformation = tgt::mat4::identity);
    ///Copy metadata from other volumehande:
    Volume(VolumeRepresentation* const volume, const VolumeBase* vh);
    Volume(VolumeRepresentation* const volume, const MetaDataContainer* mdc);
    Volume(VolumeRepresentation* const volume, const MetaDataContainer* mdc, const std::set<VolumeDerivedData*>& derivedData);

public:

    /**
     * Delete all Volume pointers and the hardware specific ones, if they have been generated.
     */
    virtual ~Volume();


    //---- meta data ----

    /**
     * Returns a container storing the meta data items
     * attached to this volume.
     */
    virtual const MetaDataContainer& getMetaDataContainer() const;

    /**
     * @overload
     */
    virtual MetaDataContainer& getMetaDataContainer();

    virtual std::vector<std::string> getMetaDataKeys() const;

    virtual const MetaDataBase* getMetaData(const std::string& key) const;

    virtual bool hasMetaData(const std::string& key) const;

    template<typename T, typename U>
    void setMetaDataValue(const std::string& key, U value);

    ///Set the MD5 hash. Should only be called by a reader.
    virtual void setHash(const std::string& hash) const;

    /// Specifies the voxel dimensions of the volume.
    virtual void setSpacing(const tgt::vec3 spacing);

    virtual void setOffset(const tgt::vec3 offset);

    virtual void setPhysicalToWorldMatrix(const tgt::mat4& transformationMatrix);

    void setModality(Modality modality);
    void setRealWorldMapping(RealWorldMapping rwm);

    void setTimestep(float timestep);


    //---- volume data and representations ----

    /// re-implemented to call VolumeBase method to prevent gcc compile errors
    template <class T>
    const T* getRepresentation() const;

    /// re-implemented to call VolumeBase method to prevent gcc compile error
    template <class T>
    void removeRepresentation();

    template <class T>
    T* getWritableRepresentation();

    /**
     * Invalidates cached values (e.g. min/max), should be called when the volume was modified.
     */
    void invalidate();

    virtual size_t getNumRepresentations() const;

    virtual const VolumeRepresentation* useConverter(const RepresentationConverterBase* converter) const;

    virtual void addRepresentation(VolumeRepresentation* rep);

    /**
     * Gives up ownership of associated representations without deleting them.
     * Calls this in order to prevent deletion of the data on destruction of the volume.
     */
    void releaseAllRepresentations();

protected:

    friend class VolumeRAM;
    friend class VolumeDiskRaw;
    friend class SliceHelper;
    friend class RawVolumeReader;

    friend class VolumeSeriesSource;
#ifdef VRN_MODULE_EXPERIMENTAL
    friend class VolumeStreamProcessor;
#endif

    /**
     * Internal getter for a volume representation with a specific index. Does not notify the VolumeMemoryManager.
     */
    virtual const VolumeRepresentation* getRepresentation(size_t i) const;

    /// Internal method for removing a representation
    virtual void removeRepresentation(size_t i);

    ///Delete all other representations.
    template<class T>
    void makeRepresentationExclusive();

    // parameter is used in destructor to prevent virtual function calls in destructor
    void deleteAllRepresentations(bool notifyObserver = true);

    /**
     * (Re)Sets the data for this volume and deletes the previous one.
     * Usually there should be no need for using this method as the volume
     * is initialized within the ctor, but some VolumeReaders need to modify
     * the read data.
     */
    void setVolumeData(VolumeRAM* const volume);

    /**
     * Reloads the volume from its origin, usually from the
     * hard disk, and regenerates the dependent hardware volumes.
     *
     * @note The Volume object as well as the dependent hardware volume objects
     *       are replaced during this operation.
     *
     * After a successful reload, volumeChanged() is called on the registered observers.
     * In case the reloading failed, the Volume's state remains unchanged.
     *
     * @return true, if the volume could be successfully reloaded.
     */
    bool reloadVolume();

    /// internal method, notifies all observers that a representation is about to be removed / deleted
    virtual void notifyRepresentationDelete(const VolumeRepresentation* rep);


    //---- member attributes ----

    mutable std::vector<VolumeRepresentation*> representations_;  ///< the mutex for the representations has to be a member within the base class, since it also has to be used there

    MetaDataContainer metaData_;

    static const std::string loggerCat_;

private:

    friend class KeyValueFactory;

    Volume(const Volume&) {} // private copy-constructor to prevent copying, use clone() instead.

    /// used in constructor to determine the OpenGL texture information from the format strings of the representation
    void determineOpenGLTypes(std::string format);

};


/////////////////////// template definitions //////////////////////////

//--- meta data ----

template<typename T, typename U>
void Volume::setMetaDataValue(const std::string& key, U value) {
    MetaDataContainer& mdc = getMetaDataContainer();
    notifyPendingDataInvalidation();
    if(mdc.hasMetaData(key)) {
        MetaDataBase* mdb = mdc.getMetaData(key);
        T* md = dynamic_cast<T*>(mdb);

        if(md)
            md->setValue(value);
        else {
            LWARNING("MetaData type mismatch! Replacing.");
            mdc.removeMetaData(key);

            T* md = new T();
            md->setValue(value);
            mdc.addMetaData(key, md);
        }
    }
    else {
        T* md = new T();
        md->setValue(value);
        mdc.addMetaData(key, md);
    }
}


//---- representation data ----

// re-implemented from VolumeBase class to prevent compile errors with gcc
template <class T>
const T* Volume::getRepresentation() const {
    return VolumeBase::getRepresentation<T>();
}

// re-implemented from VolumeBass class to prevent compile errors with gcc
template <class T>
void Volume::removeRepresentation() {
    VolumeBase::removeRepresentation<T>();
}

template <class T>
T* Volume::getWritableRepresentation() {
    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    stopRunningThreads();
    // no need to notify memory manager as it is already notified in getRepresentation
    T* rep = const_cast<T*>(getRepresentation<T>());
    if (!rep)
        return 0;
    makeRepresentationExclusive<T>();
    clearDerivedData();
    return rep;
}

template<class T>
void Volume::makeRepresentationExclusive() {

    // volumes must lock the memory manager if they use it to prevent deadlocks in multi-threading
    boost::recursive_mutex* vmmMutex = 0;
    if (VolumeMemoryManager::isInited())
        vmmMutex = VolumeMemoryManager::getRef().getMutex();
    VolumeLockGuard vmmGuard(vmmMutex);

    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    if(!VolumeBase::hasRepresentation<T>()) {
        //we would be without representations if we delete all...
        if(!getRepresentation<T>())
            return;
    }
    stopRunningThreads();

    for (auto it = representations_.begin(); it != representations_.end();) {
        VolumeRepresentation* representation = *it;

        T* test = dynamic_cast<T*>(representation);
        if (!test) {
            if (VolumeMemoryManager::isInited()) {
                if (dynamic_cast<VolumeRAM*>(representation))
                    VolumeMemoryManager::getRef().updateMainMemory();
                else if (dynamic_cast<VolumeGL*>(representation))
                    VolumeMemoryManager::getRef().updateGraphicsMemory();
            }

            delete representation;
            it = representations_.erase(it);
        }
        else
            it++;
    }
}


/*
 * Position volume centered at (0,0,0), max edge length = 1
 * WARNING: Destroys correct spacing!
 */
void VRN_CORE_API oldVolumePosition(Volume* vh);

///Center volume in world coordinates by modifying the offset.
void VRN_CORE_API centerVolume(Volume* vh);

} // namespace

#endif // VRN_VOLUME_H
