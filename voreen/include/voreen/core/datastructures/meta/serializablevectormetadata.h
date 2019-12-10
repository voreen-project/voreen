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

#ifndef VRN_SERIALIZABLEVECTORMETADATA_H
#define VRN_SERIALIZABLEVECTORMETADATA_H

#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

/**
 * The @c SerializableVectorMetaData class stores a vector of selected entities (e.g. processors).
 * In the case of processors the template parameter would be @see Processor
 *
 * @see MetaDataBase
 */

template<class T>
class SerializableVectorMetaData : public MetaDataBase {
private:
    /**
     * Help structes to distinguish between pointer and non pointer templates
     */
    template<class U>
    struct is_pointer { static const bool value = false; };

    template<class U>
    struct is_pointer<U*> { static const bool value = true; };

    template<bool>
    struct pointers{};
    /**
     * Delete and Clone functions
     */
    void deleteElements(pointers<true>) {
        //delete pointer
        for(size_t i = 0; i < values_.size(); i++)
            delete values_[i];
    }
    void deleteElements(pointers<false>) {
        //nothing to do here
    }

    MetaDataBase* cloneElements(pointers<true>) const {
        if(isOwner_) {
            std::vector<T> tmp;
            for(size_t i = 0; i < values_.size(); i++) {
                tmp.push_back(static_cast<T>(values_[i]->clone()));
            }
            return new SerializableVectorMetaData<T>(tmp,true);
        } else
            return new SerializableVectorMetaData<T>(values_,isOwner_);
     }
     MetaDataBase* cloneElements(pointers<false>) const {
         return new SerializableVectorMetaData<T>(values_,isOwner_);
     }

public:
    SerializableVectorMetaData<T>(const std::vector<T>& values = std::vector<T>(), bool isOwner = false) : values_(values), isOwner_(isOwner) {}

    virtual ~SerializableVectorMetaData() {
        if(isOwner_){
            //template magic
            deleteElements(pointers<is_pointer<T>::value>());
        }
    }

    virtual std::string getClassName() const { return "SerializableVectorMetaData<T>"; }
    virtual MetaDataBase* create() const {
        return new SerializableVectorMetaData<T>();
    }
    virtual MetaDataBase* clone() const {
        //template magic
        return cloneElements(pointers<is_pointer<T>::value>());
    }

    virtual std::string toString() const {
        return "";
    }

    virtual std::string toString(const std::string& /*component*/) const {
        return toString();
    }

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const {
        s.serialize("isOwner", isOwner_);
        s.serialize("values", values_);
    }

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s) {
        s.deserialize("isOwner",isOwner_);
        s.deserialize("values", values_);
    }

    /**
     * Sets the selected value.
     *
     * @param value The selected value
     */
    void setValues(const std::vector<T>& value) {
        if(isOwner_){
            deleteElements(pointers<is_pointer<T>::value>());
        }
        values_ = value;
    }

    /**
     * Returns the selected value.
     *
     * @return the selected value
     */
    std::vector<T> getValues() const {
        return values_;
    }

    /**
     * Adds a single selection entity
     *
     * @param value The entity
     */
    void addValue(const T& value) {
        values_.push_back(value);
    }


private:
    std::vector<T> values_;
    bool isOwner_; //< if the vector is the owner, all elements will be destroyed
};

} // namespace

#endif // VRN_SERIALIZABLEVECTORMETADATA_H
