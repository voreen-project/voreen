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

#ifndef VRN_LOCALSIMILARITYANALYSIS_H
#define VRN_LOCALSIMILARITYANALYSIS_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include "../ports/ensembledatasetport.h"

namespace voreen {


struct LocalSimilarityAnalysisInput {
    PortDataPointer<EnsembleDataset> ensemble;
    const VolumeBase* referenceVolume;
    std::unique_ptr<VolumeRAM_Float> outputVolume;
    std::string field;
    float time;
};

struct LocalSimilarityAnalysisOutput {
    std::unique_ptr<VolumeBase> volume;
};

/**
 *
 */
class VRN_CORE_API LocalSimilarityAnalysis : public AsyncComputeProcessor<LocalSimilarityAnalysisInput, LocalSimilarityAnalysisOutput>  {
public:
    LocalSimilarityAnalysis();
    virtual ~LocalSimilarityAnalysis();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "LocalSimilarityAnalysis"; }
    virtual std::string getCategory() const       { return "Plotting";                }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }

protected:

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

    void adjustPropertiesToInput();

protected:

    virtual void setDescriptions() {
        setDescription("");
    }

    EnsembleDatasetPort ensembleInport_;
    VolumePort referencePort_;
    VolumePort outport_;

    StringOptionProperty selectedField_;
    FloatProperty time_;

    static const std::string loggerCat_;
};

} // namespace

#endif
