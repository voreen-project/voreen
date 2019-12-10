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

#ifndef VRN_COLLECTIVESETTINGSPROPERTY_H
#define VRN_COLLECTIVESETTINGSPROPERTY_H

#include "voreen/core/properties/optionproperty.h"

#include <map>
#include <string>

namespace voreen {

/**
 * This property is used to define property settings modes: Each mode is mapped to
 * an arbitrary set of property settings, i.e., property/value pairs, which are applied on mode activation.
 */
class VRN_CORE_API CollectiveSettingsProperty : public StringOptionProperty {
public:
    CollectiveSettingsProperty(const std::string& id, const std::string& guiText, int invalidationLevel = Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    CollectiveSettingsProperty();
    virtual ~CollectiveSettingsProperty();

    virtual Property* create() const { return new CollectiveSettingsProperty(); }

    virtual std::string getClassName() const       { return "CollectiveSettingsProperty"; }
    virtual std::string getTypeDescription() const { return "CollectiveSettings"; }

    virtual void invalidate();

    /**
     * Adds a property settings consisting of a property/value pair to the specified mode.
     * If the mode does not exists yet, it is created.
     */
    template<class T>
    void addSetting(const std::string& mode, TemplateProperty<T>* const property, const T& value);

    /// Activates the passed mode, thereby applying the associated settings.
    void selectMode(const std::string& mode);

    virtual void select(const std::string& key);
    virtual void selectByKey(const std::string& key);

private:
    class Setting {
    public:
        virtual void applySetting() = 0;
        virtual ~Setting() {}
    };

    template<class T>
    class TemplateSetting : public voreen::CollectiveSettingsProperty::Setting {
    public:
        TemplateSetting(TemplateProperty<T>* property, const T& value);
        virtual void applySetting();
    private:
        TemplateProperty<T>* property_;
        T value_;
    };

    void applyModeSettings();

    /// maps from mode string to a vector of the corresponding property/value pairs
    std::map<std::string, std::vector<Setting*> > modeSettingsMap_;
};

// ----------------------------------------------------------------------------
// template implementations

template<class T>
voreen::CollectiveSettingsProperty::TemplateSetting<T>::TemplateSetting(TemplateProperty<T>* property, const T& value)  {
    tgtAssert(property, "null pointer passed");
    property_ = property;
    value_ = value;
}

template<class T>
void voreen::CollectiveSettingsProperty::TemplateSetting<T>::applySetting() {
    tgtAssert(property_, "null pointer");
    property_->set(value_);
}

template<class T>
void CollectiveSettingsProperty::addSetting(const std::string& mode, TemplateProperty<T>* const property, const T& value) {
    tgtAssert(mode != "", "mode must not be empty");
    tgtAssert(property, "null pointer passed");

    if (modeSettingsMap_.find(mode) == modeSettingsMap_.end()) {
        modeSettingsMap_.insert(std::make_pair(mode, std::vector<Setting*>()));
        addOption(mode, mode);
    }
    modeSettingsMap_[mode].push_back(new TemplateSetting<T>(property, value));
}

} // namespace voreen

#endif // VRN_COLLECTIVESETTINGSPROPERTY_H
