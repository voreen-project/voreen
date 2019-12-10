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

#ifndef VRN_VOLUMEMULTIPLIER_H
#define VRN_VOLUMEMULTIPLIER_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/properties/temppathproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/ports/volumeport.h"
#include "modules/hdf5/io/hdf5filevolume.h"

namespace voreen {

struct VolumeMultiplierInput {
    const VolumeBase& inputVolume;
    std::unique_ptr<HDF5FileVolume> outputVolume;

     VolumeMultiplierInput(const VolumeBase& inputVolume, std::unique_ptr<HDF5FileVolume>&& pOutputVolume)
         : inputVolume(inputVolume)
         , outputVolume(std::move(pOutputVolume))
     {
     }

     VolumeMultiplierInput(const VolumeMultiplierInput&) = delete;
     VolumeMultiplierInput(VolumeMultiplierInput&& old)
         : inputVolume(old.inputVolume)
         , outputVolume(old.outputVolume.release())
     {
     }
};
struct VolumeMultiplierOutput {
    std::string outputVolumeFilePath;
};

class VolumeMultiplier : public AsyncComputeProcessor<VolumeMultiplierInput, VolumeMultiplierOutput> {
public:
    VolumeMultiplier();
    virtual ~VolumeMultiplier();

    virtual std::string getClassName() const         { return "VolumeMultiplier";      }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual bool isEndProcessor() const       { return true; }
    virtual VoreenSerializableObject* create() const;

    virtual void setDescriptions() {
        setDescription("Processor that multiplies the size of the input volume by mirroring it in each of the coordinate axis directions. "
                "The mirroring behaviour is the same as GL_MIRRORED_REPEAT.");
    }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL; }

    virtual void initialize();
    virtual bool isReady() const;

    virtual VolumeMultiplierInput prepareComputeInput();
    virtual VolumeMultiplierOutput compute(VolumeMultiplierInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(VolumeMultiplierOutput output);

protected:
    void updateOutputSizeDisplay();

private:

    // Ports
    VolumePort inport_;
    VolumePort outport_;

    // General properties
    TempPathProperty outputVolumeFilePath_;
    IntProperty outputVolumeDeflateLevel_;

    IntVec3Property multiplicationFactor_;

    StringProperty outputSizeDisplay_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMEMULTIPLIER_H
