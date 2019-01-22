/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#ifndef VRN_PROBABILITYVOLUMECREATOR_H
#define VRN_PROBABILITYVOLUMECREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/progressproperty.h"

#include "../ports/ensembledatasetport.h"

namespace voreen {

struct ProbabilityVolumeCreatorInput {
    const EnsembleDataset& dataset;
    VolumeRAM_Float* volumeData;
    float resampleFactor;
};

struct ProbabilityVolumeCreatorOutput {
    VolumeBase* volume;
};

class VRN_CORE_API ProbabilityVolumeCreator : public AsyncComputeProcessor<ProbabilityVolumeCreatorInput, ProbabilityVolumeCreatorOutput> {
public:
    ProbabilityVolumeCreator();
    virtual ~ProbabilityVolumeCreator();

    virtual Processor* create() const;
    virtual std::string getClassName() const        { return "ProbabilityVolumeCreator"; }
    virtual std::string getCategory() const         { return "Output";                   }
    virtual CodeState getCodeState() const          { return CODE_STATE_EXPERIMENTAL;    }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:

    /// Inport for the ensemble data structure.
    EnsembleDatasetPort ensembleInport_;

    /// The ensemble data.
    VolumePort volumeOutport_;

    /// Resample factor.
    FloatProperty resampleFactor_;
};

} // namespace

#endif
