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

#ifndef VRN_VOLUMEBRICKSAVE_
#define VRN_VOLUMEBRICKSAVE_

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/volumeinfoproperty.h"

#include "modules/hdf5/io/hdf5filevolume.h"

namespace voreen {

class VolumeBrickSave : public VolumeProcessor {
public:
    VolumeBrickSave();
    virtual ~VolumeBrickSave();

    virtual std::string getClassName() const         { return "VolumeBrickSave"; }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() {
        setDescription("Processor that saves Volume to brick of a HDF5 file volume");
    }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }

    virtual void initialize();
    virtual void deinitialize();
    virtual void process();

protected:
    void saveBrick();

private:

    // Ports
    VolumePort inport_;

    // General properties
    VolumeInfoProperty volumeInfo_;
    FileDialogProperty volumeFilePath_;
    BoolProperty allowTruncateFile_;
    IntVec3Property volumeDimensions_;
    IntVec3Property brickDimensions_;
    IntVec3Property brickOffset_;
    ButtonProperty saveButton_;
    static const std::string loggerCat_;
};
} // namespace voreen

#endif // VRN_VOLUMEBRICKSAVE_
