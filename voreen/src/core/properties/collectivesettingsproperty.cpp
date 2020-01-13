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

#include "voreen/core/properties/collectivesettingsproperty.h"

namespace voreen {

CollectiveSettingsProperty::CollectiveSettingsProperty(const std::string& id, const std::string& guiText, int invalidationLevel, Property::LevelOfDetail lod)
    : StringOptionProperty(id, guiText, invalidationLevel, lod)
{ }

CollectiveSettingsProperty::CollectiveSettingsProperty()
    : StringOptionProperty()
{}

CollectiveSettingsProperty::~CollectiveSettingsProperty() {
    for (std::map<std::string, std::vector<Setting*> >::iterator it = modeSettingsMap_.begin(); it != modeSettingsMap_.end(); ++it) {
        std::vector<Setting*>& vec = it->second;
        for (size_t i=0; i<vec.size(); i++)
            delete vec.at(i);
    }
}

void CollectiveSettingsProperty::selectMode(const std::string& mode) {
    StringOptionProperty::select(mode);
    applyModeSettings();
}

void CollectiveSettingsProperty::invalidate() {
    StringOptionProperty::invalidate();
    applyModeSettings();
}

void CollectiveSettingsProperty::applyModeSettings() {
    std::string mode = get();
    if (modeSettingsMap_.find(mode) == modeSettingsMap_.end())
        return;

    std::vector<CollectiveSettingsProperty::Setting*>& settings = modeSettingsMap_[mode];
    for (size_t i=0; i<settings.size(); i++) {
        settings.at(i)->applySetting();
    }
}

void CollectiveSettingsProperty::select(const std::string& key) {
    selectMode(key);
}

void CollectiveSettingsProperty::selectByKey(const std::string& key) {
    selectMode(key);
}

}   // namespace
