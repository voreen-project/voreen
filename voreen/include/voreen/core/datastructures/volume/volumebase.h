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

#ifndef VRN_VOLUMEBASE_H
#define VRN_VOLUMEBASE_H

// representations
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumegl.h"

// derived data
#include "voreen/core/datastructures/volume/volumederiveddata.h"

// volume observer and volume url
#include "voreen/core/datastructures/volume/volumeobserver.h"
#include "voreen/core/datastructures/volume/volumeurl.h"

// volume modalities
#include "voreen/core/datastructures/volume/modality.h"

// volume memory management
#include "voreen/core/memorymanager/volumememorymanager.h"

// command queue
#include "voreen/core/utils/commandqueue.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/voreenapplication.h"

// meta data
#include "voreen/core/datastructures/meta/metadatacontainer.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

#include "voreen/core/ports/port.h"

namespace voreen {

class VolumeDerivedDataThreadBase;
class Volume;

/*
 * Helper struct to delegate VolumeObserver calls to DataInvalidationObservable
 * You probably should not be using this yourself.
 */
struct VRN_CORE_API VolumeBaseDataInvalidator : public VolumeObserver {
    VolumeBaseDataInvalidator(VolumeBase& base);

    // Methods from VolumeObserver. Used to delegate events to DataInvalidationObservable
    virtual void volumeDelete(const VolumeBase* source);
    virtual void volumeChange(const VolumeBase* source);
    virtual void volumeRepresentationDelete(const VolumeBase* source, const VolumeRepresentation* rep);
    virtual void volumeDataDelete(const VolumeBase* source);
    VolumeBase& base_;
};

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API Observable<VolumeObserver>;
#endif

/**
 * Base class for volumes. This class provides the basic interface for accessing the volume's actual data representations, meta data, and derived data.
 *
 * The VolumeDecorator class and the Volume class are derived from this class and constitute the volume class hierarchy,
 * while the actual data is stored in representations (see Volume class for more details).
 */
class VRN_CORE_API VolumeBase : public Observable<VolumeObserver>, public DataInvalidationObservable {
public:

    //---- meta data names ----
    static const std::string META_DATA_NAME_OFFSET;
    static const std::string META_DATA_NAME_SPACING;
    static const std::string META_DATA_NAME_TRANSFORMATION;
    static const std::string META_DATA_NAME_TIMESTEP;
    static const std::string META_DATA_NAME_MODALITY;
    static const std::string META_DATA_NAME_REAL_WORLD_MAPPING;

    //---- constructor and destructor ----

    VolumeBase();
    virtual ~VolumeBase();

    virtual Volume* clone() const;


    //---- observer notifications ----

    /**
     * Notifies the registered VolumeObservers about the pending
     * deletion of the Volume.
     */
    void notifyDelete();

    /**
     * Notifies the registered VolumeObservers that a reload
     * of the volume was done.
     */
    void notifyChanged();


    //---- volume origin / url ----

    /**
     * Returns the origin the volume has been loaded from,
     * usually a file path.
     */
    const VolumeURL& getOrigin() const;

    /// @overload
    VolumeURL& getOrigin();

    /**
     * Sets the origin the volume has been loaded from,
     * usually a file path.
     */
    void setOrigin(const VolumeURL& origin);


    //---- meta data ----

    virtual std::vector<std::string> getMetaDataKeys() const = 0;
    virtual const MetaDataBase* getMetaData(const std::string& key) const = 0;
    virtual bool hasMetaData(const std::string& key) const = 0;

    /*
     * Returns a meta data value for the key parameter.
     *
     * @param key The meta data key to look for.
     * @param def Default return value in case metadata could not be found.
     */
    template<typename T, typename U>
    U getMetaDataValue(const std::string& key, U def) const;

    // meta data shortcuts

    /// Returns the associated timestep of this volume.
    virtual float getTimestep() const;

    tgt::vec3 getSpacing() const;

    /**
     * The offset is specified in physical coordinates (i.e., including the spacing) and specifies the offset to the first voxel position (i.e., center of the voxel).
     */
    tgt::vec3 getOffset() const;

    virtual Modality getModality() const;

    /**
     * Used to map normalized intensity values of the data set to real-world units (e.g., Hounsfield units for CT)
     */
    RealWorldMapping getRealWorldMapping() const;


    //---- derived data ----

    /**
     * Returns the derived data item of the specified type T,
     * which must be a concrete subtype of VolumeDerivedData.
     *
     * If no derived data item of the type T exists, a new item is created
     * and stored if possible. Otherwise, 0 is returned.
     *
     * @see hasDerivedData
     */
    template<class T>
    T* getDerivedData() const;

    /**
     * Launches a thread to calculate derived data of type T.
     * Observe this volume to get notified when calculations finish (@see VolumeObserver::derivedDataThreadFinished())
     *
     * @return If a derived data of type T has already been calculated it is returned, otherwise 0.
     */
    template<class T>
    T* getDerivedDataThreaded() const;

    /**
     * Returns whether there exists a derived data item of the specified type T,
     * which must be a concrete subtype of VolumeDerivedData.
     */
    template<class T>
    T* hasDerivedData() const;

    /**
     * Adds the given data item to the derived data associated with this volume.
     * The template type T must be a concrete subtype of VolumeDerivedData.
     *
     * @note The volume takes ownership of the passed data item.
     * @note An existing item of the type T is replaced and deleted.
     */
    template<class T>
    void addDerivedData(T* data);

    void addDerivedData(VolumeDerivedData* data);

    /**
     * Removes and deletes the derived data item with the specified type T,
     * which must be a concrete subtype of VolumeDerivedData.
     * If no item with the specified type T exists, the call has no effect.
     */
    template<class T>
    void removeDerivedData();

    /**
     * Deletes all derived data items associated with this volume.
     */
    void clearDerivedData();

    template<class T>
    void derivedDataThreadFinished(VolumeDerivedDataThreadBase* ddt) const;


    //---- hashes ----

    /// Computes the MD5 hash of the raw volume data.
    virtual std::string getRawDataHash() const;
    /// Computes the MD5 hash of the meta data.
    virtual std::string getMetaDataHash() const;
    /// Concatenation of raw data and meta data hash (getRawDataHash() + "-" + getMetaDataHash())
    virtual std::string getHash() const;


    //---- data and representations

    template <class T>
    bool hasRepresentation() const;

     /*
      * Returns a representation (@see VolumeRepresentation) of type T, performing an automatic conversion between representation if necessory.
      * Cannot be used to convert between formats (e.g., from 16- to 8-bit), call using the base class (i.e., VolumeRAM).
      *
      * @return returns 0 if the representation could not be created (e.g., because of memory constraints)
      */
    template <class T>
    const T* getRepresentation() const;

    /**
     * Returns true if the volume has a representation that is NOT equal to the requested one, but can directly be converted to it.
     */
    template <class T>
    bool canConvertToRepresentation() const;

    /**
     * Tries to remove a representation (@see VolumeRepresentation) of type T.
     * If none is found, the function has no effect.
     * If one was found but it's the only representation of the volume, a warning is emitted.
     * Otherwise, the representation of type T is removed.
     */
    template <class T>
    void removeRepresentation();

    virtual size_t getNumRepresentations() const = 0;

    virtual const VolumeRepresentation* useConverter(const RepresentationConverterBase* converter) const = 0;

    virtual void addRepresentation(VolumeRepresentation* rep) = 0;

    /**
     * Get a copy of a slice in the xy-plane of this volume.
     *
     * @param sliceNumber position of the slice in z direction
     * @note sliceNumber must be < getDimensions().z
     *
     * @returns The slice as a VolumeRAM. The caller takes ownership.
     */
    VolumeRAM* getSlice(size_t sliceNumber) const;

    //---- basic volume information ----

    /// Returns the format of the volume as string (e.g., "uint8" or "Vector3(float)", @see VolumeFactory).
    virtual std::string getFormat() const;

    /// Returns the base type (e.g.,    "float" for a representation of format "Vector3(float)").
    virtual std::string getBaseType() const;

    /// Returns the OpenGL internal texture format
    virtual GLint getOpenGLInternalFormat() const;

    /// Returns the OpenGL texture format
    virtual GLenum getOpenGLFormat() const;

    /// Returns the OpenGL texture type
    virtual GLenum getOpenGLType() const;

    size_t getNumChannels() const;
    tgt::svec3 getDimensions() const;
    size_t getNumVoxels() const;
    size_t getBytesPerVoxel() const;


    //---- geometry and transformation matrices ----

    /// Returns the 8 cube vertices in physical coordinates.
    virtual std::vector<tgt::vec3> getCubeVertices() const;

    /**
     * Returns volume's bounding box as MeshGeometry.
     *
     * @param applyTransformation if true, the bounding box
     *  is transformed into world coordinates. Otherwise,
     *  the bounding box is returned in physical coordinates.
     *  @see getVoxelToWorldMatrix
     *
     * @note The mesh is internally created on each call. The bounding box includes 0.5 voxels at each border since a voxel is sampled in its center (can thus directly be used for rendering using OpenGL 3d textures).
     */
    virtual MeshGeometry getBoundingBox(bool applyTransformation = true) const;

    /// Returns the size of the cube in physical coordinates. Includes 0.5 voxels before and after the actual samples, since a voxel is sampled in the center.
    virtual tgt::vec3 getCubeSize() const;

    /// Returns the lower left front in physical coordinates. Includes 0.5 voxels before the first sample, since a voxel is sampled in its center.
    virtual tgt::vec3 getLLF() const;

    /// Returns the upper right back in physical coordinates. Includes 0.5 voxels after the last sample, since a voxel is sampled in its center.
    virtual tgt::vec3 getURB() const;

    /**
     * Returns the matrix mapping from voxel coordinates (i.e. [0; dim-1])
     * to world coordinates.
     *
     * @note The matrix is internally created on each call.
     */
    virtual tgt::mat4 getVoxelToWorldMatrix() const;

    /**
     * Returns the matrix mapping from world coordinates
     * to voxel coordinates (i.e. [0; dim-1]).
     *
     * @note The matrix is internally created on each call.
     */
    virtual tgt::mat4 getWorldToVoxelMatrix() const;

    /**
     * Returns the matrix mapping from world coordinates
     * to texture coordinates (i.e. [0.0; 1.0]).
     *
     * @note The matrix is internally created on each call.
     */
    virtual tgt::mat4 getWorldToTextureMatrix() const;

    /**
     * Returns the matrix mapping from texture coordinates (i.e. [0.0; 1.0])
     * to world coordinates.
     *
     * @note The matrix is internally created on each call.
     */
    virtual tgt::mat4 getTextureToWorldMatrix() const;

    virtual tgt::mat4 getVoxelToPhysicalMatrix() const;
    virtual tgt::mat4 getPhysicalToVoxelMatrix() const;

    virtual tgt::mat4 getPhysicalToWorldMatrix() const;
    virtual tgt::mat4 getWorldToPhysicalMatrix() const;

    virtual tgt::mat4 getTextureToPhysicalMatrix() const;
    virtual tgt::mat4 getPhysicalToTextureMatrix() const;

    virtual tgt::mat4 getTextureToVoxelMatrix() const;
    virtual tgt::mat4 getVoxelToTextureMatrix() const;


protected:

    /**
     * This class is a workaround for locking the MemoryManager mutex, which will not work using boost::lock_guard<recursive_lock>, since
     *      if (VolumeMemoryManager::isInited())
     * encapsulates it and destroys it directly afterwards as it runs out of scope. We therefore use a pointer here as a workaround.
     */
    class VolumeLockGuard {
        public:
            explicit VolumeLockGuard(boost::recursive_mutex* m = 0)
                : m_(m)
            {
                if (m_)
                    m_->lock();
            }

            virtual ~VolumeLockGuard() {
                if (m_)
                    m_->unlock();
            }

        private:
            boost::recursive_mutex* m_;
    };

    friend class VolumeDecoratorIdentity; // this is necessary to allow the VolumeDecorator to access the protected member function getRepresentation(i)

    template<class T>
        void addDerivedDataInternal(T* data) const;

    template<class T>
        void removeDerivedDataInternal() const;

    void clearFinishedThreads();
    void stopRunningThreads();

    /**
     * Internal getter for a volume representation with a specific index. Does not notify the VolumeMemoryManager.
     */
    virtual const VolumeRepresentation* getRepresentation(size_t i) const = 0;

    /// Internal method for removing a representation
    virtual void removeRepresentation(size_t i) = 0;

    VolumeURL origin_;

    // we need the mutex for the representations within the base class, even if the representations are handled by the Volume subclass!
    mutable boost::recursive_mutex representationMutex_;

    mutable std::set<VolumeDerivedData*> derivedData_; // use mutex derivedDataMutex_ to make derived data thread safe!
    mutable boost::recursive_mutex derivedDataMutex_;  // use recursive mutex to fix problems with addDerivedDataInternal which calls removeDerivedDataInternal

    mutable std::set<VolumeDerivedDataThreadBase*> derivedDataThreads_;
    mutable std::set<VolumeDerivedDataThreadBase*> derivedDataThreadsFinished_;
    mutable boost::mutex derivedDataThreadMutex_;

    // format information
    std::string format_;
    std::string baseType_;
    GLint glInternalFormat_;
    GLenum glFormat_;
    GLenum glType_;

    VolumeBaseDataInvalidator dataInvalidator_;

    static const std::string loggerCat_;
};


// template definitions -------------------------------------------------------
// ----------------------------------------------------------------------------

//---- meta data ----

template<typename T, typename U>
U VolumeBase::getMetaDataValue(const std::string& key, U def) const {
    if(hasMetaData(key)) {
        const MetaDataBase* mdb = getMetaData(key);
        const T* md = dynamic_cast<const T*>(mdb);
        if(md)
            return static_cast<U>(md->getValue());
        else
            return def;
    }
    else
        return def;
}

//---- derived data ----

template<class T>
void VolumeBase::removeDerivedData() {
    removeDerivedDataInternal<T>();
}

template<class T>
void VolumeBase::addDerivedData(T* data) {
    addDerivedDataInternal<T>(data);
}

//---- volume representations ----

template <class T>
const T* VolumeBase::getRepresentation() const {

    // volumes must lock the memory manager if they use it to prevent deadlocks in multi-threading
    boost::recursive_mutex* vmmMutex = 0;
    if (VolumeMemoryManager::isInited())
        vmmMutex = VolumeMemoryManager::getRef().getMutex();
    VolumeLockGuard vmmGuard(vmmMutex);

    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);

    // notify memory manager that this volume is in use
    if (VolumeMemoryManager::isInited())
        VolumeMemoryManager::getRef().notifyUse(const_cast<VolumeBase*>(this));

    if (getNumRepresentations() == 0) {
        LWARNING("Found no representations for this volume!" << this);
        return 0;
    }

    //Check if rep. is available:
    for (size_t i=0; i<getNumRepresentations(); i++) {
        if (const T* rep = dynamic_cast<const T*>(getRepresentation(i))) {
            return rep;
        }
    }

    // representation is not available -> since it has to be created, we will first request the memory from the VolumeMemoryManager
    if (VolumeMemoryManager::isInited()) {
        bool memoryAvailable = true;
        if (typeid(T) == typeid(VolumeRAM)) {
            // check for RAM memory
            memoryAvailable = VolumeMemoryManager::getRef().requestMainMemory(this);
            if (!memoryAvailable)
                LERROR("VolumeMemoryManager could not fulfill memory request for VolumeRAM representation");
        }
        else if (typeid(T) == typeid(VolumeGL)) {
            // check for graphics memory
            memoryAvailable = VolumeMemoryManager::getRef().requestGraphicsMemory(this);
            if (!memoryAvailable)
                LERROR("VolumeMemoryManager could not fulfill memory request for VolumeGL representation");
        }

        // if the memory manager can not allocate enough memory, the representation can not be allocated
        if (!memoryAvailable)
            return 0;
    }

    // try to convert from VolumeRAM
    ConverterFactory fac;
    if (hasRepresentation<VolumeRAM>()) {
        const VolumeRAM* volumeRam = getRepresentation<VolumeRAM>(); // it is ok to use getRepresentation<T> here, since the VolumeMemoryManager has been notified anyway
        RepresentationConverter<T>* converter = fac.findConverter<T>(volumeRam);
        if (converter) {
            const T* rep = static_cast<const T*>(useConverter(converter)); //we can static cast here because we know the converter returns T*
            if (rep)
                return rep;
        }

    }

    // try to convert from other representations
    for (size_t i=0; i<getNumRepresentations(); i++) {
        RepresentationConverter<T>* converter = fac.findConverter<T>(getRepresentation(i));
        if (converter) {
            const T* rep = static_cast<const T*>(useConverter(converter)); //we can static cast here because we know the converter returns T*
            if (rep)
                return rep;
        }
    }

    // As a last resort, try to get a RAM representation and convert using
    // that.  This of course only has a chance of working if we are not trying
    // to get a VolumeRAM itself now (for which the conversion failed earlier
    // already if we got to here).
    if (VolumeMemoryManager::isInited() && !std::is_same<VolumeRAM, T>()) {
        // Implicitly (try to) create VolumeRAM representation by recursive call.
        const VolumeRAM* volumeRam = getRepresentation<VolumeRAM>(); // it is ok to use getRepresentation<T> here, since the VolumeMemoryManager has been notified anyway
        if(volumeRam) {
            RepresentationConverter<T>* converter = fac.findConverter<T>(volumeRam);
            if (converter) {
                const T* rep = static_cast<const T*>(useConverter(converter)); //we can static cast here because we know the converter returns T*
                if (rep) {
                    return rep;
                }
            }
        } else {
            LERROR("Failed to create intermediate VolumeRAM representation");
        }
    }

    LERROR("Failed to create representation of the requested type!");
    return 0;

}

template <class T>
bool VolumeBase::hasRepresentation() const {
    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    for(size_t i=0; i<getNumRepresentations(); i++) {
        if(dynamic_cast<const T*>(getRepresentation(i)))
            return true;
    }
    return false;
}

template <class T>
void VolumeBase::removeRepresentation() {
    // volumes must lock the memory manager if they use it to prevent deadlocks in multi-threading
    boost::recursive_mutex* vmmMutex = 0;
    if (VolumeMemoryManager::isInited())
        vmmMutex = VolumeMemoryManager::getRef().getMutex();
    VolumeLockGuard vmmGuard(vmmMutex);

    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    for (size_t i = getNumRepresentations(); i > 0; i--) {
        if (dynamic_cast<const T*>(getRepresentation(i-1))) {
            // If we found a match but it's the only representation, don't delete it.
            if(getNumRepresentations() == 1) {
                LWARNING("Can't remove only representation of this volume!" << this);
                return;
            }

            // no need to notify memory manager, since removeRepresentation() already does
            removeRepresentation(i-1);
        }
    }
}

template <class T>
bool VolumeBase::canConvertToRepresentation() const {
    boost::lock_guard<boost::recursive_mutex> lock(representationMutex_);
    // try to convert from VolumeRAM (but not for VolumeRAM itself)
    ConverterFactory fac;
    if (hasRepresentation<VolumeRAM>() && (typeid(T) != typeid(VolumeRAM))) {
        // find the VolumeRAM without calling getRepresentation<VolumeRAM> (would notify the memory manager)
        const VolumeRAM* volumeRam = 0;
        for (size_t i = 0; i < getNumRepresentations(); ++i) {
            volumeRam = dynamic_cast<const VolumeRAM*>(getRepresentation(i));
            if (volumeRam)
                break;
        }
        tgtAssert(volumeRam, "Volume has VolumeRAM representation, but it cannot be found");

        RepresentationConverter<T>* converter = fac.findConverter<T>(volumeRam);
        if (converter)
            return true;
    }

    // try to convert from other representations
    for (size_t i=0; i<getNumRepresentations(); i++) {
        if (dynamic_cast<const T*>(getRepresentation(i)))
            continue;
        RepresentationConverter<T>* converter = fac.findConverter<T>(getRepresentation(i));
        if (converter)
            return true;
    }
    return false;
}

} // namespace

#include "voreen/core/datastructures/volume/volumederiveddatathread.h"

namespace voreen {

template<class T>
T* VolumeBase::getDerivedData() const {
    T* test = hasDerivedData<T>();
    if(test)
        return test;

    // Look for running threads calculating this type of derived data:
    derivedDataThreadMutex_.lock();
    for (std::set<VolumeDerivedDataThreadBase*>::iterator it=derivedDataThreads_.begin(); it!=derivedDataThreads_.end(); ++it) {
        if (typeid(**it) == typeid(VolumeDerivedDataThread<T>)) {
            VolumeDerivedDataThreadBase* ddt = *it;
            derivedDataThreadMutex_.unlock();
            if(ddt->isRunning())
                ddt->join();
            return hasDerivedData<T>();
        }
    }

    // no running thread...start a new one:
    VolumeDerivedDataThread<T>* thread = new VolumeDerivedDataThread<T>();
    derivedDataThreads_.insert(thread);
    thread->startThread(this);

    derivedDataThreadMutex_.unlock();

    thread->join(); // ...and wait for it to finish
    return hasDerivedData<T>();
}

template<class T>
T* VolumeBase::getDerivedDataThreaded() const {
    T* test = hasDerivedData<T>();
    if(test)
        return test;

    // Look for running threads calculating this type of derived data:
    derivedDataThreadMutex_.lock();
    for (std::set<VolumeDerivedDataThreadBase*>::iterator it=derivedDataThreads_.begin(); it!=derivedDataThreads_.end(); ++it) {
        if (typeid(**it) == typeid(VolumeDerivedDataThread<T>)) {
            derivedDataThreadMutex_.unlock(); // there is a thread calculating this type of derived data
            return 0;
        }
    }

    VolumeDerivedDataThread<T>* thread = new VolumeDerivedDataThread<T>();
    derivedDataThreads_.insert(thread);
    thread->startThread(this);

    derivedDataThreadMutex_.unlock();

    return 0;
}

template<class T>
void VolumeBase::derivedDataThreadFinished(VolumeDerivedDataThreadBase* ddt) const {
    VolumeDerivedData* result = ddt->getResult();
    if (result)
        addDerivedDataInternal<T>(static_cast<T*>(result));

    derivedDataThreadMutex_.lock();
    derivedDataThreads_.erase(ddt);
    derivedDataThreadMutex_.unlock();

    // enqueue notification of observers
    VoreenApplication::app()->getCommandQueue()->enqueue(this, LambdaFunctionCallback([this] {
        std::vector<VolumeObserver*> observers = Observable<VolumeObserver>::getObservers();
        for (size_t i = 0; i<observers.size(); ++i)
            observers[i]->derivedDataThreadFinished(this);
    }));

    derivedDataThreadMutex_.lock();
    derivedDataThreadsFinished_.insert(ddt);
    derivedDataThreadMutex_.unlock();
}


template<class T>
T* VolumeBase::hasDerivedData() const {
    derivedDataMutex_.lock();
    for (std::set<VolumeDerivedData*>::const_iterator it=derivedData_.begin(); it!=derivedData_.end(); ++it) {
        if (typeid(**it) == typeid(T)) {
            T* tmp = static_cast<T*>(*it);
            derivedDataMutex_.unlock();
            return tmp;
        }
    }
    derivedDataMutex_.unlock();
    return 0;
}

template<class T>
void VolumeBase::addDerivedDataInternal(T* data) const {
    if (!dynamic_cast<VolumeDerivedData*>(data)) {
        LERROR("derived data item is not of type VolumeDerivedData");
        throw std::invalid_argument("passed data item is not of type VolumeDerivedData");
    }

    derivedDataMutex_.lock();

    if (hasDerivedData<T>())
        removeDerivedDataInternal<T>();

    derivedData_.insert(static_cast<VolumeDerivedData*>(data));
    derivedDataMutex_.unlock();
}

template<class T>
void VolumeBase::removeDerivedDataInternal() const {
    if (!hasDerivedData<T>())
        return;

    derivedDataMutex_.lock();
    //T* data = getDerivedData<T>();
    for (std::set<VolumeDerivedData*>::iterator it=derivedData_.begin(); it!=derivedData_.end(); ++it) {
        if (dynamic_cast<T*>(*it)) {
            delete *it;
            derivedData_.erase(it);
            derivedDataMutex_.unlock();
            return;
        }
    }
    derivedDataMutex_.unlock();
}

} // namespace

#endif // VRN_VOLUMEBASE_H
