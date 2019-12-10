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

#ifndef VRN_TEMPPATHPROPERTY_H
#define VRN_TEMPPATHPROPERTY_H

#include "voreen/core/properties/filedialogproperty.h"

namespace voreen {

class VRN_CORE_API TempPathProperty : public FileDialogProperty {
public:

    TempPathProperty(const std::string& id, const std::string& guiText, const std::string& dialogCaption,
        const std::string& directory, const std::string& fileFilter = "", FileMode fileMode = OPEN_FILE,
        int invalidationLevel = Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    TempPathProperty();
    virtual ~TempPathProperty() {}
    virtual Property* create() const;

    virtual std::string getClassName() const { return "TempPathProperty"; }
    virtual std::string getTypeDescription() const { return "FileDialog"; }

    virtual void setUseGeneratedPath(bool use);
    virtual bool getUseGeneratedPath() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

protected:

    void generateAndUseNewTmpPath();

    bool useGeneratedPath_;
};

}   // namespace

#endif
