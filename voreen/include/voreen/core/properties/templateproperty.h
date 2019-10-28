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

#ifndef VRN_TEMPLATEPROPERTY_H
#define VRN_TEMPLATEPROPERTY_H

#include "voreen/core/properties/property.h"
#include "voreen/core/properties/templatepropertycondition.h"
#include "voreen/core/datastructures/callback/callbackmanager.h"

namespace voreen {

template<class T>
class PropertyMutator;

/**
 * Template for a property that stores a single value.
 */
template<class T>
class TemplateProperty : public Property {
public:
    TemplateProperty(const std::string& id, const std::string& guiText,
                     T value, int = Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    TemplateProperty();

    //TemplateProperty(const TemplateProperty*);

    virtual ~TemplateProperty();

    virtual void set(const T& value);

    virtual const T& get() const { return value_; }

    // Allow access value access to a mutator of this property
    friend T& PropertyMutator<T>::getValue() const;

    PropertyMutator<T> getMutator();

    void setDefaultValue(const T& value);

    const T& getDefault() const { return defaultValue_; }

    virtual void reset();

    Condition* addValidation(const Condition& condition) {
        validations_.push_back(condition.clone());
        return validations_.back();
    }

    TemplateProperty<T>& verifies(const Condition& condition) {
        addValidation(condition);
        return *this;
    }

    /**
     * Returns whether the given value would be valid for this property.
     * The property is not modified.
     *
     * @param value the value to check
     * @param errorMsg contains the error message, if the validation failed
     */
    bool isValidValue(const T& value, std::string& errorMsg);

protected:
    /**
     * Runs validate() on all Validations. If any of them is not met() a Condition::ValidationFailed
     * exception is thrown.
     */
    void validate(const T& value, std::string& errorMsg, bool restore = true);

    T value_;
    T defaultValue_;
    std::vector<Condition*> validations_;

//private:
//    std::string getTypename() const;
};

/**
 * Class that allows mutable access to the value of a template property,
 * but ensures* that the property will be updated after mutation, i.e.
 * when the mutator is destroyed. This is useful for properties with
 * non scalar values.
 *
 * *: The caller shall not store and not access a reference obtained
 * via the mutator after the mutator has been destroyed.
 */
template<class T>
class PropertyMutator {
public:
    PropertyMutator(TemplateProperty<T>* owner);
    PropertyMutator(PropertyMutator&& other);

    //PropertyMutator(const PropertyMutator& other) = delete;
    //PropertyMutator<T>& operator=(const PropertyMutator& other) = delete;

    ~PropertyMutator();
    T* operator->() const;
    T& operator*() const;
    T& getValue() const;

private:

    // Disallow copying
    PropertyMutator(const PropertyMutator& other);
    PropertyMutator<T>& operator=(const PropertyMutator& other);

    TemplateProperty<T>* owner_;
};

// TemplateProperty implementation ---------------------------------------------------------------

class TransFuncBase;

template<class T>
TemplateProperty<T>::TemplateProperty(const std::string& id, const std::string& guiText,
                                      T value, int invalidationLevel, Property::LevelOfDetail lod)
  : Property(id, guiText, invalidationLevel, lod)
  , value_(value)
  , defaultValue_(value)
{}

template<class T>
TemplateProperty<T>::TemplateProperty() :
    Property()
{}

template<>
TemplateProperty<TransFuncBase*>::TemplateProperty();

template<class T>
TemplateProperty<T>::~TemplateProperty() {
    size_t i;

    // delete validations
    for (i = 0; i < validations_.size(); ++i)
        delete validations_.at(i);

    validations_.clear();
}

template<class T>
void TemplateProperty<T>::set(const T& value) {
    if (value_ != value) {
        std::string errorMsg;
        validate(value, errorMsg, false);

        if (!errorMsg.empty())
            LWARNINGC("voreen.TemplateProperty", errorMsg);
        if (value_ != value)
            return;

        invalidate();  // issues invalidateOwner and updateWidgets
    }
}

template<class T>
PropertyMutator<T> TemplateProperty<T>::getMutator() {
   return std::move(PropertyMutator<T>(this));
}

template<class T>
void TemplateProperty<T>::setDefaultValue(const T& value) {
    defaultValue_ = value;
}

template<class T>
void TemplateProperty<T>::reset() {
    set(defaultValue_);
}

template<class T>
void TemplateProperty<T>::validate(const T& value, std::string& errorMsg, bool restore) {
    // save value
    T temp(value_);
    value_ = value;
    bool valid = true;
    for (size_t j = 0; j < validations_.size() && valid; ++j)
        valid &= validations_[j]->validate(errorMsg);

    if (!valid || restore) // restore if new value valid or requested
        value_ = temp;
}

template<class T>
bool TemplateProperty<T>::isValidValue(const T& value, std::string& errorMsg) {
    // save value
    T temp(value_);
    value_ = value;
    bool valid = true;
    for (size_t j = 0; j < validations_.size() && valid; ++j)
        valid &= validations_[j]->validate(errorMsg);
    value_ = temp;
    return valid;
}

// PropertyMutator implementation ---------------------------------------------------------------

template<class T>
PropertyMutator<T>::PropertyMutator(TemplateProperty<T>* owner)
    : owner_(owner)
{
}
template<class T>
PropertyMutator<T>::PropertyMutator(PropertyMutator&& other)
    : owner_(other.owner_)
{
    other.owner_ = nullptr;
}

template<class T>
PropertyMutator<T>::~PropertyMutator() {
    if(owner_) {
        owner_->invalidate();
    }
}

template<class T>
T* PropertyMutator<T>::operator->() const {
    return &getValue();
}

template<class T>
T& PropertyMutator<T>::operator*() const {
    return getValue();
}

template<class T>
T& PropertyMutator<T>::getValue() const {
    return owner_->value_;
}

} // namespace voreen

#endif // VRN_TEMPLATEPROPERTY_H
