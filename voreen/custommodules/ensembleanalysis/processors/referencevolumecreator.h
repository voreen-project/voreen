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

#ifndef VRN_REFERENCEVOLUMECREATOR_H
#define VRN_REFERENCEVOLUMECREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/vectorproperty.h"

#include "../ports/ensembledatasetport.h"

namespace voreen {
    
struct ReferenceVolumeCreatorInput {
    PortDataPointer<EnsembleDataset> ensemble;
    std::unique_ptr<VolumeRAM> outputVolume;
    std::string field;
    float time;
    std::string referenceMethod;
    size_t referenceRun;
};

struct ReferenceVolumeCreatorOutput {
    std::unique_ptr<VolumeBase> volume;
};

/**
 *
 */
class VRN_CORE_API ReferenceVolumeCreator : public AsyncComputeProcessor<ReferenceVolumeCreatorInput, ReferenceVolumeCreatorOutput>  {
public:
    ReferenceVolumeCreator();
    virtual ~ReferenceVolumeCreator();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "ReferenceVolumeCreator"; }
    virtual std::string getCategory() const       { return "Ensemble Processing";    }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;  }

protected:

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

    void adjustPropertiesToInput();

protected:

    virtual void setDescriptions() {
        setDescription("Creates a reference volume from the input ensemble for a selected time step and field.");
        referenceMethod_.setDescription("Use the reference method to determine how the reference volume is created.");
    }

    EnsembleDatasetPort inport_;
    VolumePort outport_;

    StringOptionProperty selectedField_;
    FloatProperty time_;

    StringOptionProperty referenceMethod_;
    StringOptionProperty referenceRun_;

    IntVec3Property outputDimensions_;

    static const std::string loggerCat_;
};

} // namespace

#endif
