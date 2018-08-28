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

#ifndef VRN_TRANSFERFUNCPROPERTYGENERIC_H
#define VRN_TRANSFERFUNCPROPERTYGENERIC_H

#include "voreen/core/properties/transfunc/transfuncpropertybase.h"

namespace voreen {

/**
 * Generic property for all transfer functions.
 * Used to get rid of dynamic casts.
 */
template<typename T>
class TransFuncPropertyGeneric : public TransFuncPropertyBase {
public:

    /**
     * Constructor
     *
     * @param ident identifier that is used in serialization
     * @param guiText text that is shown in the gui
     * @param invalidationLevel The owner is invalidated with this InvalidationLevel upon change.
     * @param lod level of detail in the gui representation
     */
    TransFuncPropertyGeneric(const std::string& ident, const std::string& guiText, int invalidationLevel = Processor::INVALID_RESULT,
                          Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    TransFuncPropertyGeneric();
    virtual ~TransFuncPropertyGeneric();

    //-------------------------------------------------
    //  Property Functions
    //-------------------------------------------------


    /**
     * Sets the stored transfer function to the given one.
     *
     * @note The TransFuncProperty takes ownership of the passed
     *  object. Therefore, the caller must not delete it.
     *
     * @param tf transfer function the property is set to
     */
    void set(T* tf);

    /**
     * Returns the current transfer function.
     */
    virtual T* get() const;

    /**
     * Creates an initial transfer function, if has not already been created by deserialization.
     *
     * @see Property::initialize
     */
    void initialize();

    /**
     * Deletes the stored transfer function.
     *
     * @see Property::deinitialize
     */
    void deinitialize();

    /** @see Property::serialize */
    virtual void serialize(Serializer& s) const;
    /** @see Property::deserialize */
    virtual void deserialize(Deserializer& s);

    //-------------------------------------------------
    //  Member
    //-------------------------------------------------
protected:
    T* transFunc_;  ///< pointer to the current function.

    //-------------------------------------------------
    //  Hide base pointer of super class
    //-------------------------------------------------
private:
    TransFuncPropertyBase::baseFunction_;
    virtual void set(TransFuncBase* tf) {TransFuncPropertyBase::set(tf);}
};



//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
template<typename T>
TransFuncPropertyGeneric<T>::TransFuncPropertyGeneric(const std::string& ident, const std::string& guiText, int invalidationLevel, Property::LevelOfDetail lod)
    : TransFuncPropertyBase(ident, guiText, invalidationLevel, lod)
    , transFunc_(0)
{}

template<typename T>
TransFuncPropertyGeneric<T>::TransFuncPropertyGeneric()
    : TransFuncPropertyBase("", "", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , transFunc_(0)
{}

template<typename T>
TransFuncPropertyGeneric<T>::~TransFuncPropertyGeneric() {
    //delete done in base class
    transFunc_ = 0;
}



template<typename T>
void TransFuncPropertyGeneric<T>::serialize(Serializer& s) const {
    TransFuncPropertyBase::serialize(s);

    s.serialize("TransferFunction", transFunc_);
}

template<typename T>
void TransFuncPropertyGeneric<T>::deserialize(Deserializer& s) {
    TransFuncPropertyBase::deserialize(s);

    T* tf = 0;
    s.deserialize("TransferFunction", tf);
    set(tf);
}

template<typename T>
void TransFuncPropertyGeneric<T>::initialize() {
    TransFuncPropertyBase::initialize();

    // create initial transfer function, if it has not been created during deserialization
    if (!transFunc_) {
        set(new T());
        LGL_ERROR;
    }
}

template<typename T>
void TransFuncPropertyGeneric<T>::deinitialize() {
    if (transFunc_) {
        delete transFunc_;
        transFunc_ = 0;
        baseFunction_ = 0;
        LGL_ERROR;
    }

    TransFuncPropertyBase::deinitialize();
}

template<typename T>
void TransFuncPropertyGeneric<T>::set(T* tf) {
    // tf object already assigned
    if (tf == transFunc_) {
        return;
    }

    // assign new object, but store previous one for deletion
    T* oldValue = transFunc_;
    transFunc_ = tf;

    //notify changes
    TransFuncPropertyBase::set(tf);

    //delete old function (deletes base function too)
    delete oldValue;
}

template<typename T>
T* TransFuncPropertyGeneric<T>::get() const {
    return transFunc_;
}



} // namespace voreen

#endif // VRN_TRANSFERFUNCPROPERTYGENERIC_H
