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

#ifndef VRN_VOLUMESAVE_H
#define VRN_VOLUMESAVE_H

#include <string>
#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/volumeinfoproperty.h"

namespace voreen {

class Volume;
class VolumeSerializer;
class VolumeSerializerPopulator;

class VRN_CORE_API VolumeSave : public VolumeProcessor {
public:
    VolumeSave();
    virtual ~VolumeSave();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "VolumeSave";      }
    virtual std::string getCategory() const   { return "Output";          }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE; }
    virtual bool isEndProcessor() const       { return true;              }

    virtual bool usesExpensiveComputation() const { return true; }

    virtual void saveVolume();

protected:
    virtual void setDescriptions() {
        setDescription("Saves a volume to disk.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();
    virtual void invalidate(int inv = 1);

private:
    std::vector<std::string> constructFilters() const;

    VolumePort inport_;

    FileDialogProperty filename_;
    ButtonProperty saveButton_;
    BoolProperty saveOnPathChange_;
    BoolProperty continousSave_;
    VolumeInfoProperty volumeInfo_;

    bool saveVolume_;

    VolumeSerializerPopulator* volSerializerPopulator_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
