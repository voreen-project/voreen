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

#include "voreen/core/properties/temppathproperty.h"

#include "voreen/core/voreenapplication.h"

namespace voreen {

TempPathProperty::TempPathProperty(const std::string& id, const std::string& guiText,
    const std::string& dialogCaption, const std::string& directory,
    const std::string& fileFilter, FileDialogProperty::FileMode fileMode,
    int invalidationLevel, Property::LevelOfDetail lod)
    : FileDialogProperty(id, guiText, dialogCaption, directory, fileFilter, fileMode, invalidationLevel, lod, VoreenFileWatchListener::ALWAYS_OFF)
{
    // Initialize auto generated path.
    setUseGeneratedPath(true);
}

TempPathProperty::TempPathProperty()
    : FileDialogProperty()
{}

Property* TempPathProperty::create() const {
    return new TempPathProperty();
}

void TempPathProperty::setUseGeneratedPath(bool use) {
    useGeneratedPath_ = use;
    if (useGeneratedPath_ && VoreenApplication::app()) {
        generateAndUseNewTmpPath();
    }
    updateWidgets();
}

bool TempPathProperty::getUseGeneratedPath() const {
    return useGeneratedPath_;
}

void TempPathProperty::serialize(Serializer& s) const {
    FileDialogProperty::serialize(s);
    s.serialize("useGeneratedPath", useGeneratedPath_);
}

void TempPathProperty::deserialize(Deserializer& s) {
    FileDialogProperty::deserialize(s);
    s.deserialize("useGeneratedPath", useGeneratedPath_);
    if(useGeneratedPath_) {
        generateAndUseNewTmpPath();
    }
}

void TempPathProperty::generateAndUseNewTmpPath() {
    tgtAssert(VoreenApplication::app(), "No voreen application");
    set(VoreenApplication::app()->getUniqueTmpFilePath(fileFilter_));
}

}   // namespace
