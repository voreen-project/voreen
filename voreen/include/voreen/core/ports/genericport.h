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

#ifndef VRN_GENERICPORT_H
#define VRN_GENERICPORT_H

#include "voreen/core/ports/port.h"
#include "voreen/core/datastructures/imagesequence.h"
#include "voreen/core/datastructures/cloneable.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/datastructures/callback/callbackmanager.h"

#include <sstream>

namespace tgt {
    class Texture;
}

namespace voreen {


/*
 * A smart pointer that may either be a direct pointer to port data or a copy of the data in the port.
 *
 * @see getThreadSafeData()
 */
template<typename T>
class PortDataPointer {
public:
    PortDataPointer(const T* pointer, bool owned)
        : pointer_(pointer)
        , owned_(owned)
    {
    }
    PortDataPointer(PortDataPointer&& other)
        : pointer_(other.pointer_)
        , owned_(other.owned_)
    {
        other.pointer_ = nullptr;
    }
    ~PortDataPointer() {
        if(owned_) {
            delete pointer_;
        }
    }
    PortDataPointer& operator=(PortDataPointer&& other) {
        if(this != &other) {
            // Destruct the current object, but keep the memory.
            this->~PortDataPointer();
            // Call the move constructor on the memory region of the current object.
            new(this) PortDataPointer(std::move(other));
        }

        return *this;
    }
    const T& operator*() const {
        return *pointer_;
    }
    const T* operator->() const {
        return pointer_;
    }
    const T* get() const {
        return pointer_;
    }
    operator const T*() const {
        return pointer_;
    }
    operator bool() const {
        return pointer_ != nullptr;
    }

private:

    // Delete copy and copy assignment.
    PortDataPointer(const PortDataPointer&);
    PortDataPointer& operator=(const PortDataPointer&);

    const T* pointer_;
    bool owned_;
};

/**
 * @brief Template port class to store points to type T.
 *
 * The data is always stored in the outport, inports fetch data from connected outports.
 */
template<typename T>
class GenericPort : public Port {
public:
    GenericPort(PortDirection direction, const std::string& id, const std::string& guiName = "", bool allowMultipleConnections = false,
                         Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT);
    virtual ~GenericPort();

    /**
     * Returns whether GenericPort::getData() returns 0 or not,
     * therefore indicating if there is any data at the port.
     */
    virtual bool hasData() const;

    /**
     * Set data stored in this port. Can only be called on outports.
     * @param takeOwnership If true, the data will be deleted by the port.
     */
    virtual void setData(const T* data, bool takeOwnership = true);

    /// Return the data stored in this port (if this is an outport) or the data of the first connected outport (if this is an inport).
    virtual const T* getData() const;

    /// Return the data stored in this port (if this is an outport) or the data of all the connected outports (if this is an inport).
    virtual std::vector<const T*> getAllData() const;

    /// Return whether or not the data is owned and will thus be deleted by the port (see setData).
    bool ownsData() const;

    /**
     * Get a pointer to the data of the port that is safe to use in a thread IF the thread is cancelled on invalidation of the port data.
     * For processing data independently of the evaluation of the network consider AsyncComputeProcessor or refer to it as a reference
     * implementation on safe handling of port data.
     *
     * @see DataInvalidationObservable
     * @see AsyncComputeProcessor
     */
    PortDataPointer<T> getThreadSafeData() const;

    /**
     * Gets a vector from pointers as described in getThreadSafeData, hence returning a thread safe data representation
     * for each contained data.
     *
     * @see getThreadSafeData
     * @see DataInvalidationObservable
     * @see AsyncComputeProcessor
     */
    std::vector<PortDataPointer<T>> getThreadSafeAllData() const;

    /**
     * Register a callback
     * @param callback Is called if the data of the port changes and data exists
     */
    void onNewData(const Callback& callback);


    std::vector<const GenericPort<T>* > getConnected() const;

    /**
     * Returns true, if the port is connected and
     * contains a data object.
     */
    virtual bool isReady() const;

    virtual void invalidatePort();

    virtual void clear();

    virtual bool isDataInvalidationObservable() const;
    virtual bool isDataThreadSafe() const;
    virtual void addDataInvalidationObserver(const DataInvalidationObserver* observer) const;
    virtual void removeDataInvalidationObserver(const DataInvalidationObserver* observer) const;

protected:
    const T* portData_;
    bool ownsData_;
    CallbackManager onNewDataCallbacks_;
};



#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API GenericPort<VolumeList>;
#endif

class VRN_CORE_API VolumeListPort : public GenericPort<VolumeList> {

public:

    VolumeListPort(PortDirection direction, const std::string& id, const std::string& guiName = "", bool allowMultipleConnections = false, Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT)
    : GenericPort<VolumeList>(direction, id, guiName, allowMultipleConnections, invalidationLevel) {}

    virtual std::string getClassName() const {return "VolumeListPort";}
    virtual Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const {return new VolumeListPort(direction,id,guiName);}
    virtual tgt::col3 getColorHint() const {
        return tgt::col3(255, 0, 255);
    }
    virtual std::string getContentDescription() const {
        std::stringstream strstr;
        strstr << Port::getContentDescription();
        if(hasData())
            strstr << std::endl << "Size of List: " << getData()->size();
        return strstr.str();
    }
    virtual std::string getContentDescriptionHTML() const {
        std::stringstream strstr;
        strstr << Port::getContentDescriptionHTML();
        if(hasData())
            strstr << "<br>" << "Size of List: " << getData()->size();
        return strstr.str();
    }
};

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API GenericPort<ImageSequence>;
#endif

class VRN_CORE_API ImageSequencePort : public GenericPort<ImageSequence> {

public:

    ImageSequencePort(PortDirection direction, const std::string& id, const std::string& guiName = "", bool allowMultipleConnections = false, Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT)
    : GenericPort<ImageSequence>(direction, id, guiName, allowMultipleConnections, invalidationLevel) {}

    virtual Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const {return new ImageSequencePort(direction,id,guiName);}
    virtual std::string getClassName() const { return "ImageSequencePort"; }

    virtual std::string getContentDescription() const {
        std::stringstream strstr;
        strstr << Port::getContentDescription();
        if(hasData())
            strstr << std::endl << "Size of List: " << getData()->size();
        return strstr.str();
    }

    virtual std::string getContentDescriptionHTML() const {
        std::stringstream strstr;
        strstr << Port::getContentDescriptionHTML();
        if(hasData())
            strstr << "<br>" << "Size of List: " << getData()->size();
        return strstr.str();
    }

};


// ---------------------------------------- implementation ----------------------------------------

template <typename T>
GenericPort<T>::GenericPort(PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections,
                     Processor::InvalidationLevel invalidationLevel)
    : Port(direction, id, guiName, allowMultipleConnections, invalidationLevel)
    , portData_(0)
    , ownsData_(false)
{}

template <typename T>
GenericPort<T>::~GenericPort() {
    notifyDataWillChange();

    if (ownsData_) {
        delete portData_;
    }
    portData_ = nullptr;

    notifyDataHasChanged();
}

template <typename T>
void GenericPort<T>::setData(const T* data, bool takeOwnership) {
    tgtAssert(isOutport(), "called setData on inport!");

    if(data != portData_) {
        notifyDataWillChange();

        //delete previous data
        if (ownsData_)
            delete portData_;

        portData_ = data;
        notifyDataHasChanged();
    }

    ownsData_ = takeOwnership;

    invalidatePort();
}

template <typename T>
void GenericPort<T>::invalidatePort(){
    Port::invalidatePort();
    if(hasData()){
        onNewDataCallbacks_.execute();
    }

}

template <typename T>
const T* GenericPort<T>::getData() const {
    if (isOutport())
        return portData_;
    else {
        for (size_t i = 0; i < connectedPorts_.size(); ++i) {
            if (!connectedPorts_[i]->isOutport())
                continue;

            GenericPort<T>* p = static_cast< GenericPort<T>* >(connectedPorts_[i]);
            if (p)
                return p->getData();
        }
    }
    return 0;
}

template <typename T>
bool GenericPort<T>::hasData() const {
    return (getData() != nullptr);
}

// Give out data pointer for DataInvalidationObservable types.
template<typename T>
static PortDataPointer<T> getThreadSafeDataInternal(const T* data, typename std::enable_if<std::is_base_of<DataInvalidationObservable, T>::value>::type* = 0) {
    return PortDataPointer<T>(data, false);
}

// Give out data pointer for Cloneable types.
template<typename T>
static PortDataPointer<T> getThreadSafeDataInternal(const T* data, typename std::enable_if<!std::is_base_of<DataInvalidationObservable, T>::value && std::is_base_of<Cloneable<T>, T>::value>::type* = 0) {
    return PortDataPointer<T>(data ? data->clone().release() : nullptr, true);
}

// Stop compilation for types that are neither DataInvalidationObservable nor Cloneable.
template<typename T>
struct __DelayedInstatiationFalse : std::false_type
{ };
template<typename T>
static PortDataPointer<T> getThreadSafeDataInternal(const T* data, typename std::enable_if<!std::is_base_of<DataInvalidationObservable, T>::value && !std::is_base_of<Cloneable<T>, T>::value>::type* = 0) {
#ifdef _MSC_VER // FIXME: enable opt-out of unused static functions in MSVC at compile time
    tgtAssert(false, "Threadsafe generic port data must implement Cloneable or DataInvalidationObservable!");
    return PortDataPointer<T>(nullptr, false);
#else
    static_assert(__DelayedInstatiationFalse<T>::value, "Threadsafe generic port data must implement Cloneable or DataInvalidationObservable!");
#endif
}

template <typename T>
PortDataPointer<T> GenericPort<T>::getThreadSafeData() const {
    return getThreadSafeDataInternal(getData());
}

template <typename T>
std::vector<const T*> GenericPort<T>::getAllData() const {
    std::vector<const T*> allData;

    if (isOutport())
        allData.push_back(portData_);
    else {
        for (size_t i = 0; i < connectedPorts_.size(); ++i) {
            if (!connectedPorts_[i]->isOutport())
                continue;
            GenericPort<T>* p = static_cast<GenericPort<T>*>(connectedPorts_[i]);
            if(p->getData())
                allData.push_back(p->getData());
        }
    }

    return allData;
}
template <typename T>
bool GenericPort<T>::ownsData() const {
    return ownsData_;
}

template <typename T>
std::vector<PortDataPointer<T>> GenericPort<T>::getThreadSafeAllData() const {
    std::vector<PortDataPointer<T>> allDataThreadSafe;
    for (const T* data : getAllData()) {
        allDataThreadSafe.push_back(std::move(getThreadSafeDataInternal(data)));
    }
    return allDataThreadSafe;
}

template <typename T>
std::vector<const GenericPort<T>*> GenericPort<T>::getConnected() const {
    std::vector<const GenericPort<T>*> ports;
    for (size_t i = 0; i < connectedPorts_.size(); ++i) {
        GenericPort<T>* p = static_cast<GenericPort<T>*>(connectedPorts_[i]);

        ports.push_back(p);
    }
    return ports;
}

template <typename T>
bool GenericPort<T>::isReady() const {
    if (isOutport())
        return isConnected();
    else
        return (!getConnected().empty() && hasData() && areConditionsMet());
}

template <typename T>
void GenericPort<T>::clear() {
    if (isOutport()) {
        if(portData_) {
            setData(nullptr, false);
        }
    }
    else
        LERROR("clear() called on inport");
}

template<typename T>
void GenericPort<T>::onNewData( const Callback& callback ) {
    onNewDataCallbacks_.registerCallback(callback);
}

template<typename T>
bool GenericPort<T>::isDataThreadSafe() const {
    return std::is_base_of<DataInvalidationObservable, T>::value || std::is_base_of<Cloneable<T>, T>::value;
}

template<typename T>
bool GenericPort<T>::isDataInvalidationObservable() const {
    return std::is_base_of<DataInvalidationObservable, T>::value;
}

// Give out data pointer for polymorphic types
template<typename T>
static const DataInvalidationObservable* getInvalidationObservableDataInternal(const GenericPort<T>& p, typename std::enable_if<std::is_base_of<DataInvalidationObservable, T>::value && std::is_polymorphic<T>::value>::type* = 0) {
    return dynamic_cast<const DataInvalidationObservable*>(p.getData());
}

// Give out data pointer for non-polymorphic types.
template<typename T>
static const DataInvalidationObservable* getInvalidationObservableDataInternal(const GenericPort<T>& p, typename std::enable_if<std::is_base_of<DataInvalidationObservable, T>::value && !std::is_polymorphic<T>::value>::type* = 0) {
    return static_cast<const DataInvalidationObservable*>(p.getData());
}

// In case T does not derive from DataInvalidationsObservable, the function must not be called either way and Port::add/removeDataInvalidationObserver will fail with an assertion.
template<typename T>
static const DataInvalidationObservable* getInvalidationObservableDataInternal(const GenericPort<T>& p, typename std::enable_if<!std::is_base_of<DataInvalidationObservable, T>::value>::type* = 0) {
    tgtAssert(false, "DataInvalidationObservers can only be added/removed if generic port data implements DataInvalidationObservable!");
    return nullptr;
}

template<typename T>
void GenericPort<T>::addDataInvalidationObserver(const DataInvalidationObserver* observer) const {
    Port::addDataInvalidationObserver(observer);
    getInvalidationObservableDataInternal(*this)->addObserver(observer); // Safe, since Port::addDataInvalidationObserver checks validity.
}

template<typename T>
void GenericPort<T>::removeDataInvalidationObserver(const DataInvalidationObserver* observer) const {
    Port::removeDataInvalidationObserver(observer);
    getInvalidationObservableDataInternal(*this)->removeObserver(observer); // Safe, since Port::removeDataInvalidationObserver checks validity.
}

} // namespace

#endif // VRN_GENERICPORT_H
