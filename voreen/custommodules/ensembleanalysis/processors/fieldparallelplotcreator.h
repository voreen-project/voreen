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

#ifndef VRN_FIELDPARALLELPLOTCREATOR_H
#define VRN_FIELDPARALLELPLOTCREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"

#include "../datastructures/fieldplotdata.h"
#include "../ports/ensembledatasetport.h"
#include "../ports/fieldplotdataport.h"
#include "voreen/core/ports/volumeport.h"


namespace voreen {

struct FieldParallelPlotCreatorInput {
    const EnsembleDataset& dataset;
    std::unique_ptr<FieldPlotData> outputPlot;
    std::vector<tgt::vec3> seedPoints;
};

struct FieldParallelPlotCreatorOutput {
    std::unique_ptr<FieldPlotData> plotData;
};

/**
 *
 */
class VRN_CORE_API FieldParallelPlotCreator : public AsyncComputeProcessor<FieldParallelPlotCreatorInput, FieldParallelPlotCreatorOutput> {
public:
    static const std::string META_DATA_HASH;


    FieldParallelPlotCreator();
    virtual ~FieldParallelPlotCreator();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "FieldParallelPlotCreator"; }
    virtual std::string getCategory() const       { return "Plotting";                 }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;    }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:

    virtual bool isReady() const;
    virtual void adjustPropertiesToInput();
    virtual std::vector<std::reference_wrapper<Port>> getCriticalPorts();

    virtual void setDescriptions() {
        setDescription("");
    }

    EnsembleDatasetPort inport_;
    VolumePort seedMask_;
    FieldPlotDataPort outport_;

    IntProperty numSeedPoints_;
    IntProperty seedTime_;
    IntProperty verticalResolution_;
    IntProperty horizontalResolutionPerTimeUnit_;

    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_FIELDPARALLELPLOTCREATOR_H
