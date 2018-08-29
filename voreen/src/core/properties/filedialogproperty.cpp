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

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/condition.h"
#include "voreen/core/utils/voreenfilepathhelper.h"

#include "tgt/filesystem.h"

namespace voreen {

const std::string FileDialogProperty::loggerCat_("voreen.FileDialogProperty");

FileDialogProperty::FileDialogProperty(const std::string& id, const std::string& guiText,
                               const std::string& dialogCaption, const std::string& directory,
                               const std::string& fileFilter, FileDialogProperty::FileMode fileMode,
                               int invalidationLevel, Property::LevelOfDetail lod,
                               VoreenFileWatchListener::WatchMode watchMode)
    : StringProperty(id, guiText, "", invalidationLevel, lod)
    , VoreenFileWatchListener(watchMode)
    , dialogCaption_(dialogCaption)
    , directory_(directory)
    , fileMode_(fileMode)
    , fileFilter_(fileFilter)
    , recursiveSet_(false)
{
    if (fileMode == DIRECTORY)
        set(directory_);
}

FileDialogProperty::FileDialogProperty()
    : StringProperty()
    , VoreenFileWatchListener(VoreenFileWatchListener::OPTIONAL_OFF)
    , fileMode_(OPEN_FILE)
    , recursiveSet_(false)
{}

Property* FileDialogProperty::create() const {
    return new FileDialogProperty();
}

void FileDialogProperty::fileActionCallback() {
    invalidate();
}

const std::string& FileDialogProperty::getDialogCaption() const {
    return dialogCaption_;
}

void FileDialogProperty::setDialogCaption(const std::string& caption) {
    dialogCaption_ = caption;
}

const std::string& FileDialogProperty::getDirectory() const {
    return directory_;
}

void FileDialogProperty::setDirectory(const std::string& dir) {
    directory_ = dir;
}

const std::string& FileDialogProperty::getFileFilter() const {
    return fileFilter_;
}

void FileDialogProperty::setFileFilter(const std::string& filter) {
    fileFilter_ = filter;
}

FileDialogProperty::FileMode FileDialogProperty::getFileMode() const {
    return fileMode_;
}

void FileDialogProperty::setFileMode(FileMode mode) {
    fileMode_ = mode;
}

void FileDialogProperty::set(const std::string& value) {
    if (get() == value)
        return;

    // Determine if set was called recursively.
    // We only want to modify the file watch on the first recursion level.
    bool recursion = recursiveSet_;
    if (!recursiveSet_) {
        recursiveSet_ = true;
        removeWatch(get());
    }

    StringProperty::set(value);

    if (!recursion) {
        bool success = addWatch(get());
        if (!success) {
            LWARNING("Parent directory of " << tgt::FileSystem::fileName(value) << " does not exist. Resetting path.");
            setFileWatchEnabled(false);
            StringProperty::set("");
        }
        recursiveSet_ = false;
    }
}

void FileDialogProperty::serialize(Serializer& s) const {
    Property::serialize(s);
    VoreenFileWatchListener::serialize(s);

    s.serialize("paths",VoreenFilePathHelper(value_));
}

void FileDialogProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);
    VoreenFileWatchListener::deserialize(s);

    VoreenFilePathHelper tmp;
    try {
        s.deserialize("paths", tmp);
    }
    catch (SerializationNoSuchDataException) {
        s.removeLastError();
        LINFO("trying old deserialization");
        //old deserialization
        try {
            bool noPathSet;
            s.deserialize("noPathSet", noPathSet);
            if(noPathSet) {
                set("");
                return;
            }
        }
        catch (SerializationNoSuchDataException) {
            s.removeLastError();
        }
        std::string value;
        s.deserialize("value", value);

        // convert path relative to the document's path to an absolute one
        if (!s.getDocumentPath().empty() && !tgt::FileSystem::isAbsolutePath(value))
            value = tgt::FileSystem::absolutePath(tgt::FileSystem::dirName(s.getDocumentPath()) + "/" + value);

        try {
            set(value);
        }
        catch (Condition::ValidationFailed& e) {
            s.addError(e);
        }
        return;
    }
    // if no valid path has been found
    catch (VoreenException& e) {
        LERROR(e.what());
    }

    try {
        set(tmp.getPath());
    }
    catch (Condition::ValidationFailed& e) {
        s.addError(e);
    }
}

}   // namespace
