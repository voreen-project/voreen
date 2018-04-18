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

#ifndef VRN_FILEDIALOGPROPERTY_H
#define VRN_FILEDIALOGPROPERTY_H

#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/utils/voreenfilewatcher.h"

namespace voreen {

class VRN_CORE_API FileDialogProperty : public StringProperty, public VoreenFileWatchListener {
public:

    enum FileMode {
        OPEN_FILE = 0,
        DIRECTORY = 1,
        SAVE_FILE = 2
    };

    FileDialogProperty(const std::string& id, const std::string& guiText, const std::string& dialogCaption,
                   const std::string& directory, const std::string& fileFilter = "", FileMode fileMode = OPEN_FILE,
                   int invalidationLevel=Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT,
                   VoreenFileWatchListener::WatchMode watchMode = VoreenFileWatchListener::OPTIONAL_OFF);
    FileDialogProperty();
    virtual ~FileDialogProperty() {}
    virtual Property* create() const;

    virtual std::string getClassName() const       { return "FileDialogProperty"; }
    virtual std::string getTypeDescription() const { return "FileDialog"; }

    const std::string& getDialogCaption() const;
    void setDialogCaption(const std::string& caption);

    const std::string& getDirectory() const;
    void setDirectory(const std::string& dir);

    const std::string& getFileFilter() const;
    void setFileFilter(const std::string& filter);

    FileMode getFileMode() const;
    void setFileMode(FileMode mode);

    virtual void set(const std::string& value);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    /**
     * Inherited from VoreenFileWatchListener.
     * This implementation invalidates the property.
     */
    virtual void fileActionCallback() override;

protected:
    std::string dialogCaption_;
    std::string directory_;
    FileMode fileMode_;
    std::string fileFilter_;
    bool recursiveSet_;
};

}   // namespace

#endif
