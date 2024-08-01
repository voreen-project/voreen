/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEOCTREELEVELEXTRACTOR_H
#define VRN_VOLUMEOCTREELEVELEXTRACTOR_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"

namespace voreen {

class VRN_CORE_API VolumeOctreeLevelExtractor : public VolumeProcessor {

public:
    VolumeOctreeLevelExtractor();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Volume Processing"; }
    virtual std::string getClassName() const  { return "VolumeOctreeLevelExtractor";     }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING;   }

protected:

    virtual void setDescriptions() {
        setDescription("Extracts one level of a volume octree and outputs it as a VolumeRAM. Does not modify the input volume in any way during the process.");
        enableProcessing_.setDescription("If not set, the processor will not perform any volume extraction, even if the corresponding button is pressed.");
        mode_.setDescription("Selects if the resolution of the extracted level is determined by the maximum edge length in voxels or by the octree level.");
        //channel_.setDescription("Sets the channel for which the volume is extracted.");
        maxEdgeLength_.setDescription("Sets the maximum edge length for the extracted volume which can be used to determine the octree level. For the maximum edge length 2^k, the size of the largest dimension of the volume will be in the closed interval [ 2^(k-1) + 1, 2^k], depending on the actual input volume size.");
        level_.setDescription("Selectd the octree level on which the brick buffers are used to create the output volume.");
        extractButton_.setDescription("Creates the output volume (unless the processor is disabled)");
        automaticallyComputeOnChange_.setDescription("If enabled, a new output volume is computed every time the input volume or settings of the processor change.");
    }

    virtual void initialize();
    //virtual void deinitialize();

    virtual void process();

    /// adjust the property options (e.g. max octree level) to the input data set
    virtual void adjustPropertiesToInput();

    /// adjust which property is set to read-only according to the selected mode, also enables or disables the button
    virtual void adjustPropertyActivationState();

    /// synchronize the octree level and max edge length properties
    virtual void synchronizePropertiesOnChange();

    /// called when the button has been pressed, sets a flag that the volume extraction has to be computed
    virtual void buttonPressed();

private:

    VolumePort inport_;
    VolumePort outport_;

    BoolProperty enableProcessing_;
    StringOptionProperty mode_;
    //IntProperty channel_;
    IntOptionProperty maxEdgeLength_;
    IntProperty level_;
    ButtonProperty extractButton_;
    BoolProperty automaticallyComputeOnChange_;

    bool buttonPressed_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif // VRN_VOLUMEOCTREELEVELEXTRACTOR_H
