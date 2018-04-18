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

#ifndef VRN_VOLUMEBRICKSOURCE_H
#define VRN_VOLUMEBRICKSOURCE_H

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/volumeinfoproperty.h"

#include "modules/hdf5/io/hdf5filevolume.h"
#include "voreen/core/datastructures/volume/volumedisk.h"

namespace voreen {

class VolumeBrickSource : public VolumeProcessor {
public:
    VolumeBrickSource();
    virtual ~VolumeBrickSource();

    virtual std::string getClassName() const         { return "VolumeBrickSource"; }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() {
        setDescription("Processor that loads a brick of a HDF5 input volume");
    }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }

    virtual void initialize();
    virtual void deinitialize();
    virtual void process();

protected:
    void tryOpenVolume();
    void loadBrick();
    void updateBrickDimensions();
    void updateCurrentBrickOffset();
    bool isValidNumBricks(const tgt::ivec3&) const;
    template<int dim>
    bool isValidNumBricksInDim(const tgt::ivec3&) const;
    template<int dim>
    tgt::ivec3 makeValidNumBricks(const tgt::ivec3&) const;
    void adaptToChangedNumBricks();

    void setOutput(const VolumeBase* volume);

    tgt::svec3 getCurrentBrickDimensions() const;

private:

    // Ports
    VolumePort outport_;

    // General properties
    FileDialogProperty volumeFilePath_;
    VolumeInfoProperty volumeInfo_;
    IntVec3Property numBricks_;
    IntVec3Property volumeDimensions_;
    IntVec3Property brickDimensions_;
    IntVec3Property currentBrickOffset_;
    IntProperty brickToLoad_;
    ButtonProperty loadButton_;
    VolumeInfoProperty brickInfo_;

    std::unique_ptr<VolumeBase> currentVolume_;
    VolumeDisk* currentVolumeDisk_;

    static const std::string loggerCat_;
};
} // namespace voreen

#endif // VRN_VOLUMEBRICKSOURCE_H
